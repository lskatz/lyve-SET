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
  GetOptions($settings,qw(help reference=s tempdir=s altFreq=s coverage=i));

  my $bam=$ARGV[0] or die "ERROR: need bam\n".usage();
  $samplename=basename($bam,qw(.sorted.bam .bam));
  my $reference=$$settings{reference} or die "ERROR: need --reference\n".usage();
  $$settings{tempdir}||="tmp";
  $$settings{coverage}||=10;
  $$settings{altFreq}||=0.75;

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
  system("samtools mpileup -f '$reference' '$bam' 1>$pileup");
  die "ERROR: problem with\n  samtools mpileup -f '$reference' '$bam'" if $?;
  return $pileup;
}

sub varscan{
  my($pileup,$settings)=@_;
  die "ERROR: the pileup is a zero-byte file\n  $pileup" if(-s $pileup < 1);
  my $vcf="$$settings{tempdir}/".fileparse($pileup).".tmp.vcf.gz";
  system("varscan.sh mpileup2cns $pileup --min-coverage $$settings{coverage} --min-var-freq $$settings{altFreq} --output-vcf 1 |\
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

  $vcf->add_header_line({key=>'FORMAT', ID=>'AF', Number=>'A', Type=>'Float', Description=>"Allele Frequency"});
  $vcf->add_header_line({key=>'FORMAT', ID=>'RF', Number=>'A', Type=>'Float', Description=>"Reference Frequency"});
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

    # Put in allele frequency for the one accepted ALT
    # Alternate calls are in the AD field for varscan
    $$x{gtypes}{$samplename}{AD}||=0;
    $$x{gtypes}{$samplename}{RD}||=0;
    $$x{gtypes}{$samplename}{AF}||=0;
    $$x{gtypes}{$samplename}{RF}||=0;
    if($$x{gtypes}{$samplename}{AD} > 0){
      $$x{gtypes}{$samplename}{AF}=$$x{gtypes}{$samplename}{AD}/$$x{gtypes}{$samplename}{DP};
      $$x{gtypes}{$samplename}{RF}=$$x{gtypes}{$samplename}{RD}/$$x{gtypes}{$samplename}{DP};
    }

    #$vcf->add_format_field($x,"AC");
    #$$x{gtypes}{$samplename}{AC}||=$$x{gtypes}{$samplename}{AD};

    # Mask low-quality sites
    if( $$depth{$posId} < $$settings{coverage} ||             # low coverage
        (length($$x{ALT}[0]) > 1 || length($$x{REF}) > 1)     # indel
      ){
      $$x{ALT}[0]="N";
    }

    # Mask unsupported base calls
    #if( 
    #    ($is_snp && $$x{gtypes}{$samplename}{AF} < $$settings{altFreq})
    #    ||
    #    (!$is_snp && $$x{gtypes}{$samplename}{RF} < $$settings{altFreq})
    #  ){
    #  $$x{ALT}[0]="N";
    #}

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


sub usage{
  local $0=fileparse $0;
  "$0: find SNPs using varscan
  Usage: $0 file.bam > file.vcf
  -t tmp/ A temporary directory to store files
  -c 10   Min coverage
  -a 0.75 Min consensus agreement for a SNP
  "
}
