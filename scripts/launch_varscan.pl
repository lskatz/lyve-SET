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

use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/@fastaExt @fastqExt @bamExt logmsg/;
use Vcf;

# I want to make $samplename global in this script,
# which is a rare decision for me. Its value will be
# changed later.
my $samplename="SampleNameUnknown";

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
  $samplename=basename($bam,@bamExt);                    # remove bam extension
  $samplename=~s/\-$refname$//;                          # remove reference genome name
  $samplename=basename($samplename,@fastqExt,@fastaExt); # remove fast[aq] extension
  $samplename=basename($samplename,@fastaExt);           # remove fasta extension

  my $vcf=varscan2($bam,$reference,$settings);

  # This script pipes to stdout, so go for it.
  system("cat $vcf");
  die "ERROR: missing tmp vcf file $vcf, or I could not read it" if $?;

  system("rm -f $vcf");

  return 0;
}

sub varscanWorker{
  my($bam,$reference,$regionQueue,$settings)=@_;

  my $TID=threads->tid;

  # Figure out mpileup options
  my $mpileupxopts="";
  $mpileupxopts.="-Q 0 -B "; # assume that all reads have been properly filtered at this point and that mappings are good
  $mpileupxopts.="--positions $$settings{region} " if($$settings{region});

  my $tmpdir=tempdir("$$settings{tempdir}/$samplename/varscanWorker.XXXXXX");

  my $vcfCounter=0;
  while(defined(my $region=$regionQueue->dequeue)){
    my $pileup="$tmpdir/TID$TID.$vcfCounter.pileup";
    my $vcf="$tmpdir/TID$TID.$vcfCounter.vcf";
    # Avoid disk I/O problems by sleeping
    # a bit of time, depending on the number of threads
    sleep ($TID % $$settings{numcpus});

    # Run mpileup on this one region
    my $mpileupCommand="samtools mpileup -f $reference $mpileupxopts --region '$region' $bam > $pileup";
    system($mpileupCommand);
    die "ERROR running mpileup on region $region:\n  $mpileupCommand" if $?;

    # Estimate the size of the pileup. If it's zero, then
    # there will be no SNPs in it and should be skipped.
    if( -s $pileup){
      logmsg "There is a pileup in region $region and so I will run varscan";
    } else {
      logmsg "No pileup was formed in region $region and so there will be no SNPs. Skipping.";
      next;
    }

    # Finally run varscan on the pileup. Bgzip and index the
    # output file. Later it will be combined with bcftools concat.
    my $varscanCommand="varscan.sh mpileup2cns $pileup --min-coverage $$settings{coverage} --min-var-freq $$settings{altFreq} --output-vcf 1 --min-avg-qual 0 > $vcf.tmp.vcf";
    system($varscanCommand);
    die "ERROR with varscan!\n  $varscanCommand" if $?;

    system("bgzip $vcf.tmp.vcf && tabix $vcf.tmp.vcf.gz");
    die if $?;
    
    # Fix the VCF and rename the sample
    logmsg "Fixing the VCF into a new file $vcf";
    system("set_fixVcf.pl --numcpus 1 --min_coverage $$settings{coverage} --min_alt_frac $$settings{altFreq} --fail-samples --fail-sites --DP4 2 --rename-sample $samplename $vcf.tmp.vcf.gz > $vcf");
    die if $?;

    # Compress and index
    system("bgzip $vcf"); die if $?;
    system("tabix $vcf.gz"); die if $?;

    # Save the VCF
    system("mv -v $vcf.gz $vcf.gz.tbi $$settings{tempdir}/$samplename/ >&2");
    die if $?;

    $vcfCounter++;
  }

  # Some cleanup
  system("rm -rf $tmpdir");
}

sub varscan2{
  my($bam,$reference,$settings)=@_;
  # I like to use the shorthand $b for basename, but it is a reserved
  # variable name in sort{} and so I need to use 'local' instead of 'my'
  local $b=fileparse($bam);
  my $pileup="$$settings{tempdir}/$b.mpileup";
  return $pileup if(-e $pileup && -s $pileup > 0);

  # Make sure this whole operation is self-contained
  # and that the temporary dir exists for this sample.
  my $tmpdir="$$settings{tempdir}/$samplename";
  system("rm -rf $tmpdir");
  mkdir($tmpdir);
  # Make sure that no other pileup file gets in the way.
  system("rm -vf $$settings{tempdir}/$b*.mpileup $$settings{tempdir}/$b*.vcf.gz* >&2");

  # Kick off multithreading
  my $regionQueue=Thread::Queue->new();
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&varscanWorker,$bam,$reference,$regionQueue,$settings);
  }

  # Figure out any regions
  my @regions=`makeRegions.pl --numcpus $$settings{numcpus} --numchunks $$settings{numcpus} $bam`;
  chomp(@regions);
  # Add the regions and some thread-terminators to the queue
  $regionQueue->enqueue(@regions);
  $regionQueue->enqueue(undef) for(@thr);
  # When each thread encounters an undef, it will terminate.
  # When the threads are terminated and joined, we
  # can continue past this line.
  for(@thr){
    $_->join;
    die "ERROR with thread ".$_->tid if($_->error());
  }

  # Sort all subVCFs by their starting position since each one is 
  # already internally sorted.
  my @subVcf=glob("$tmpdir/*.vcf.gz");
  @subVcf=sort{
    my($chromA,$posA)=firstVcfPos($a);
    my($chromB,$posB)=firstVcfPos($b);
    return $chromA cmp $chromB if($chromA ne $chromB);
    return $posA <=> $posB;
  } @subVcf;
  my $subVcf=join(" ",@subVcf);

  # Join all pieces of the sample's vcf. All pieces
  # have been sorted at this point.
  system("bcftools concat $subVcf > $$settings{tempdir}/$b.merged.vcf");
  die if $?;

  # Make sure everything is cleaned up whenever the script ends
  system("rm -rf $tmpdir");

  return "$$settings{tempdir}/$b.merged.vcf";
}

# Return the first CHROM and POS in a given vcf.gz.
sub firstVcfPos{
  my($vcf,$settings)=@_;
  my($CHROM,$POS) = ("",0);

  open(VCF,"zcat $vcf | ") or die "ERROR: could not open $vcf: $!";
  while(<VCF>){
    next if(/^#/);
    ($CHROM,$POS)=split(/\t/);
    last if($CHROM);
  }
  close VCF;
  return($CHROM,$POS);
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
