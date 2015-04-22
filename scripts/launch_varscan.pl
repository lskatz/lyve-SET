#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Creates a vcf using varscan

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename fileparse/;
use Bio::Perl;

use FindBin;
use lib "$FindBin::RealBin/../lib/vcftools_0.1.12b/perl";
use Vcf;

my $samplename="Sample1"; # to be changed later
$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());
 

sub main{
  my $settings={};
  GetOptions($settings,qw(help reference=s tempdir=s altFreq=s coverage=i region=s)) or die $!;

  my $bam=$ARGV[0] or die "ERROR: need bam\n".usage();
  $samplename=basename($bam,qw(.sorted.bam .bam));
  my $reference=$$settings{reference} or die "ERROR: need --reference\n".usage();
  $$settings{tempdir}||="tmp";
  $$settings{coverage}||=10;
  $$settings{altFreq}||=0.75;
  $$settings{region}||=0;

  # Check to see if varscan is installed correctly
  `varscan.sh >/dev/null`;
  die "ERROR: varscan.sh gave an error" if $?;

  my $pileup=mpileup($bam,$reference,$settings);
  my $vcf=varscan($pileup,$settings);
  backfillVcfValues($vcf,$bam,$settings);

  # remove temporary files
  unlink($_) for($pileup,$vcf);

  return 0;
}

sub mpileup{
  my($bam,$reference,$settings)=@_;
  my $pileup="$$settings{tempdir}/".fileparse($bam).".mpileup";
  return $pileup if(-e $pileup && -s $pileup > 0);
  logmsg "Creating a pileup $pileup";
  
  my $xopts="";
  $xopts.="-Q 0 -B "; # assume that all reads have been properly filtered at this point and that mappings are good
  $xopts.="--positions $$settings{region} " if($$settings{region});
  my $command="samtools mpileup -f '$reference' $xopts '$bam' 1>$pileup";
  logmsg "Running mpileup:\n  $command";
  system($command);
  die "ERROR: problem with mpileup" if $?;
  return $pileup;
}

sub varscan{
  my($pileup,$settings)=@_;
  die "ERROR: the pileup is a zero-byte file\n  $pileup" if(-s $pileup < 1);
  my $vcf="$$settings{tempdir}/".fileparse($pileup).".tmp.vcf.gz";
  system("varscan.sh mpileup2cns $pileup --min-coverage $$settings{coverage} --min-var-freq $$settings{altFreq} --output-vcf 1 --min-avg-qual 0 |\
    perl -lane 's/Sample1/\Q$samplename\E/; print;' |\
    bgzip -c > $vcf
  ");
  die "ERROR: problem with either varscan.sh, bgzip, or perl one-liner for regex" if $?;

  return $vcf;
}

sub backfillVcfValues{
  my($vcfFile,$bam,$settings)=@_;
  logmsg "Backfilling values in $vcfFile and printing to stdout";
  my $vcf=Vcf->new(file=>$vcfFile);

  # Reasons a SNP doesn't pass
  my $fail_lowCoverage="DP$$settings{coverage}";
  my $fail_lowRefFreq ="RF$$settings{altFreq}";
  my $fail_lowAltFreq ="AF$$settings{altFreq}";
  my $fail_indel      ="isIndel";
  
  # Add new headers
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowCoverage, Description=>"Depth is less than $$settings{coverage}, the user-set coverage threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowRefFreq, Description=>"Reference variant consensus is less than $$settings{altFreq}, the user-set threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowAltFreq, Description=>"Allele variant consensus is less than $$settings{altFreq}, the user-set threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_indel, Description=>"Indels are not used for analysis in Lyve-SET"});

  # Genotype fields
  $vcf->add_header_line({key=>'FORMAT', ID=>'AF', Number=>'A', Type=>'Float', Description=>"Allele Frequency"});
  $vcf->add_header_line({key=>'FORMAT', ID=>'RF', Number=>'A', Type=>'Float', Description=>"Reference Frequency"});
  $vcf->add_header_line({key=>'FORMAT', ID=>'ADP', Number=>'A', Type=>'Integer', Description=>"Depth of bases with Phred score >= 15"});

  # Remove unwanted headers
  $vcf->remove_header_line(key=>'INFO',ID=>'ADP');

  # Figure out the coverage depth
  my $depth=covDepth($bam,$settings);

  # Done with headers; parse them
  $vcf->parse_header();

  # start printing
  print $vcf->format_header();
  while(my $x=$vcf->next_data_hash()){
    my $posId=$$x{CHROM}.':'.$$x{POS};
    my $pos=$$x{POS};
    $$x{gtypes}{$samplename}{DP}=$$depth{$posId} || 0;

    # Only one call is allowed for ALT
    $$x{ALT}=[$$x{ALT}[0]];
    $$x{FILTER}=[$$x{FILTER}[0]];

    # Put in exactly what the alternate call is
    $$x{ALT}[0]=$$x{REF} if($$x{ALT}[0] eq '.');
    my $is_snp=0;
    if($$x{ALT}[0] ne $$x{REF}){
      $is_snp=1;
    }

    # some info fields belong in format fields (deleted in the next step from INFO)
    for(qw(ADP)){
      $$x{gtypes}{$samplename}{$_}=$$x{INFO}{$_};
    }
    # I just don't think some fields even belong here
    for(qw(ADP HET NC HOM WT)){
      delete($$x{INFO}{$_});
    }

    # Put in allele frequency for the one accepted ALT
    # Alternate calls are in the AD field for varscan
    $$x{gtypes}{$samplename}{AD}||=0;
    $$x{gtypes}{$samplename}{RD}||=0;
    $$x{gtypes}{$samplename}{AF}||=0;
    $$x{gtypes}{$samplename}{RF}||=0;
    if($$x{gtypes}{$samplename}{ADP} > 0){
      $$x{gtypes}{$samplename}{AF}=$$x{gtypes}{$samplename}{AD}/$$x{gtypes}{$samplename}{ADP};
      $$x{gtypes}{$samplename}{RF}=$$x{gtypes}{$samplename}{RD}/$$x{gtypes}{$samplename}{ADP};
    }
    # round the frequencies
    $_=sprintf("%0.3f",$_) for($$x{gtypes}{$samplename}{AF},$$x{gtypes}{$samplename}{RF});

    ###########
    # MASKING
    # Mask low coverage
    if( $$depth{$posId} < $$settings{coverage} ){
      vcf_markAmbiguous($x,$fail_lowCoverage,$settings);
    }

    # Mask indels
    if(length($$x{ALT}[0]) > 1 || length($$x{REF}) > 1){
      vcf_markAmbiguous($x,$fail_indel,$settings);
    }

    # Mask low frequency
    if($is_snp && $$x{gtypes}{$samplename}{AF} < $$settings{altFreq}){
      vcf_markAmbiguous($x,$fail_lowRefFreq,$settings);
    }
    if(!$is_snp && $$x{gtypes}{$samplename}{RF} < $$settings{altFreq}){
      vcf_markAmbiguous($x,$fail_lowAltFreq,$settings);
    }
    # END MASKING
    ###############

    print $vcf->format_line($x);
  }

  return 1;
}

sub covDepth{
  my($bam,$settings)=@_;
  my $depthFile="$bam.depth";
  my %depth;
  die if $?;
  open(IN,"samtools depth '$bam' | ") or die "Could not open $bam in samtools:$!";
  while(<IN>){
    my($rseq,$pos,$depth)=split(/\t/);
    chomp($depth);
    $depth{$rseq.":".$pos}=$depth;
  }
  close IN;
  return \%depth;
}

# VCF line manipulations
sub vcf_markAmbiguous{
  my($x,$reason,$settings)=@_;
  $$x{ALT}[1]="N";
  $$x{gtypes}{$samplename}{GT}="2/2";
  @{$$x{FILTER}}=() if($$x{FILTER}[0] eq "PASS");
  push(@{$$x{FILTER}},$reason);
}


sub usage{
  local $0=fileparse $0;
  "$0: find SNPs using varscan
  Usage: $0 file.bam --reference ref.fasta > file.vcf
  --tempdir  tmp/      A temporary directory to store files
  --coverage 10        Min coverage
  --altFreq  0.75      Min consensus agreement for a SNP
  --region   file.bed  File of positions to read (default with no bed file: read all positions)
  "
}
