#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Creates a vcf using varscan

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename fileparse/;
use File::Temp qw/tempdir/;
use Bio::Perl;
use Bio::FeatureIO;

use FindBin;
use lib "$FindBin::RealBin/../lib/vcftools_0.1.12b/perl";
use Vcf;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/@fastaExt @fastqExt @bamExt/;

my $samplename="Sample1"; # to be changed later
$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());
 

sub main{
  my $settings={};
  GetOptions($settings,qw(help reference=s tempdir=s altFreq=s coverage=i region=s exclude=s numcpus=i)) or die $!;

  my $bam=$ARGV[0] or die "ERROR: need bam\n".usage();
  die "ERROR: only one bam allowed right now" if (@ARGV > 1);
  my $reference=$$settings{reference} or die "ERROR: need --reference\n".usage();
  my $refname=basename($reference,@fastaExt);

  $$settings{tempdir} ||=tempdir("set_varscan.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{coverage}||=10;
  $$settings{altFreq} ||=0.75;
  $$settings{region}  ||="";
  $$settings{exclude} ||="";
  $$settings{numcpus} ||=1;

  ## Parameter error checking
  # Not sure if I should include this error check or not
  # die "ERROR: Cannot use --region and --exclude together\n".usage() if($$settings{exclude} && $$settings{region});
  die "ERROR: altFreq must be between 0 and 1, inclusive" if($$settings{altFreq} < 0 || $$settings{altFreq} > 1);
  die "ERROR: coverage must be >=0" if($$settings{coverage} < 0);

  # Check to see if varscan is installed correctly
  `varscan.sh >/dev/null`;
  die "ERROR: varscan.sh gave an error" if $?;

  # Transform the sample name, e.g.,
  # NC001416.fasta.wgsim.fastq.gz-reference.sorted.bam => NC001416.fasta.wgsim
  # sample1.fastq.gz-reference.sorted.bam => sample1
  $samplename=basename($bam,@bamExt);
  $samplename=~s/\-$refname$//;
  $samplename=basename($samplename,@fastqExt,@fastaExt);
  $samplename=basename($samplename,@fastaExt);

  # To make varscan work, first do mpileup.
  # Then, it reads from mpileup.
  # Lyve-SET needs to alter some of the VCF in the third step, backfillVcfValues().
  my $pileup=mpileup($bam,$reference,$settings);
  my $vcf=varscan($pileup,$settings);
  backfillVcfValues($vcf,$bam,$settings);

  # remove temporary files
  unlink($_) for($pileup,$vcf);

  return 0;
}

sub mpileup{
  my($bam,$reference,$settings)=@_;
  my $b=fileparse($bam);
  my $pileup="$$settings{tempdir}/$b.mpileup";
  return $pileup if(-e $pileup && -s $pileup > 0);
  logmsg "Creating a pileup $pileup";

  # Figure out any defined seqnames
  my @seqname;
  open(SAMTOOLS,"samtools view -H '$bam' | ") or die "ERROR: could not open $bam with samtools: $!";
  while(<SAMTOOLS>){
    next if(!/^\@SQ/);
    chomp;
    my @F=split /\t/;
    for(@F){
      my($key,$value)=split /:/;
      push(@seqname,$value) if($key eq 'SN');
    }
  }
  close SAMTOOLS;
  my $regions=join("\n",@seqname)."\n";

  # Figure out mpileup options
  my $xopts="";
  $xopts.="-Q 0 -B "; # assume that all reads have been properly filtered at this point and that mappings are good
  $xopts.="--positions $$settings{region} " if($$settings{region});

  # Multithread a pileup so that each region gets piped into a file.
  # Then, these individual files can be combined at a later step.
  system("rm -fv $$settings{tempdir}/$b.*.mpileup >&2");
  my $command="echo \"$regions\" | xargs -P $$settings{numcpus} -n 1 -I {} sh -c 'echo \"MPileup on {}\" >&2; samtools mpileup -f $reference $xopts --region \"{}\" $bam > $$settings{tempdir}/$b.\$\$.mpileup' ";
  logmsg "Running mpileup:\n  $command";
  system($command);
  die "ERROR with xargs and samtools mpileup" if $?;

  # Concatenate the output into a single file.
  # Individual regions are already sorted, thanks to the way mpileup works.
  logmsg "Sorting mpileup results into a combined file";
  system("cat $$settings{tempdir}/$b.*.mpileup > $pileup && rm -f $$settings{tempdir}/$b.*.mpileup");
  die "ERROR with sorting mpileup results and then deleting the intermediate files" if $?;

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

  my $excludeSites=readBed($$settings{exclude},$settings) if($$settings{exclude});

  logmsg "Backfilling values in $vcfFile and printing to stdout";
  my $vcf=Vcf->new(file=>$vcfFile);

  # Reasons a SNP doesn't pass
  my $fail_lowCoverage="DP$$settings{coverage}";
  my $fail_lowRefFreq ="RF$$settings{altFreq}";
  my $fail_lowAltFreq ="AF$$settings{altFreq}";
  my $fail_indel      ="isIndel";
  my $fail_masked     ="masked";

  # How was the bam generated?
  my $mapper="unknown";
  my $CL="unknown";
  open(SAMTOOLS,"samtools view -H '$bam' | ") or die "ERROR: could not open $bam with samtools: $!";
  while(<SAMTOOLS>){
    next if(!/^\@PG/);
    chomp;
    my @F=split /\t/;
    for(@F){
      my($key,$value)=split /:/;
      $mapper=$value if($key eq "ID");
      $CL=$value if($key eq "CL");
    }
  }
  
  # Add new headers
  $vcf->add_header_line({key=>'reference',value=>$$settings{reference}});
  $vcf->add_header_line({key=>'mapper',value=>$mapper});
  $vcf->add_header_line({key=>'mapperCL',value=>$CL});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowCoverage, Description=>"Depth is less than $$settings{coverage}, the user-set coverage threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowRefFreq, Description=>"Reference variant consensus is less than $$settings{altFreq}, the user-set threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowAltFreq, Description=>"Allele variant consensus is less than $$settings{altFreq}, the user-set threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_indel, Description=>"Indels are not used for analysis in Lyve-SET"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_masked, Description=>"This site was masked using a bed file or other means"});

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

    # It's a SNP if it isn't the same as REF and isn't a dot
    my $is_snp=0;
    if($$x{ALT}[0] ne $$x{REF} && $$x{ALT}[0] ne '.'){
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
    # Mask sites that are outright masked
    if($$excludeSites{$$x{CHROM}}{$pos}){
      vcf_markAmbiguous($x,$fail_masked,$settings);
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
  my $ucRef=uc($$x{REF}); # uppercase reference, to make it easier for str comparison

  # Figure out the alt allele (N) and the correct genotype integer
  my $gtInt=0; # for the GT, e.g., 1/1 or 2/2
  # Mark if it is the reference allele
  if($ucRef eq $$x{ALT}[0] || $ucRef eq '.'){
    #shift(@{ $$x{ALT} }); # unnecessary to mark the ref base in ALT
    $gtInt=0;
  } else {
    $gtInt=1;
  }
  #push(@{ $$x{ALT} }, "N");
  $$x{ALT}=["N"];
  $$x{gtypes}{$samplename}{GT}="$gtInt/$gtInt";

  # empty the filter and then put on the fail reason
  @{$$x{FILTER}}=() if($$x{FILTER}[0] eq "PASS");
  push(@{$$x{FILTER}},$reason);
  #$$x{FILTER}=[$reason];
}

# Get a hash of seqname->{pos} from a bed file
sub readBed{
  my($bed,$settings)=@_;

  my %bed;
  my $bedin=Bio::FeatureIO->new(-format=>"bed",-file=>$bed);
  while(my $feat=$bedin->next_feature){
    my $seqname=$feat->seq_id;
    my $start=$feat->start;
    my $end=$feat->end;
    for my $pos($start..$end){
      $bed{$seqname}{$pos}=1;
    }
  }
  $bedin->close;
  return \%bed;
}

sub usage{
  local $0=fileparse $0;
  "$0: find SNPs using varscan
  Usage: $0 file.bam --reference ref.fasta > file.vcf
  --numcpus  1         How many cpus to use
  --tempdir  tmp/      A temporary directory to store files
  --coverage 10        Min coverage
  --altFreq  0.75      Min consensus agreement for a SNP
  --region   file.bed  File of positions to include (default with no bed file: read all positions)
  --exclude  file.bed  File of positions to mask    (conflicts with --region; default: don't exclude)
  "
}
