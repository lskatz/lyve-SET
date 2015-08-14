#!/usr/bin/env perl

## Simple version of the BayesHammer wrapper script.
## Assumes user provides only paired-end Illumina reads as input.
## Paired reads are provided as output, Unpaired reads discarded.
## Requires SPAdes 3.0 or higher and cg_pipeline/scripts/run_assembly_shuffleReads.pl


use strict;
use warnings;
use Data::Dumper;
use Cwd;
use Getopt::Long;
use Bio::Perl;
use File::Basename qw/basename fileparse/;
use File::Temp qw/tempdir/;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

use LyveSET qw/logmsg/;

local $0=basename($0);
exit (main());

sub main{
  my $settings={};
  
  GetOptions($settings,qw(help tempdir=s numcpus=s)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("XXXXXX",TMPDIR=>1,CLEANUP=>1);
  my @infile=@ARGV;
  die usage() if($$settings{help} || !@infile);

  logmsg "Temporary directory is $$settings{tempdir}";

  for my $infile(@infile){
    errorCorrect($infile,$settings);
  }
  
  return 0;
}

sub errorCorrect{
  my($infile,$settings)=@_;

  # Is it paired end?
  my $spadesInParam="--12 $infile";
  my $isPE=`run_assembly_isFastqPE.pl $infile`+0;
  if($isPE == 0){
    $spadesInParam="-s $infile";
  }

  # Run spades
  #logmsg "DEBUG - skipping spades";
  system("spades.py $spadesInParam -o $$settings{tempdir} --only-error-correction --disable-gzip-output -t $$settings{numcpus} >&2");
  die "ERROR with spades.py" if $?;

  # Find the correct reads to shuffle. One set of reads is
  # the set of singletons and should be discarded.
  # TODO: worry about an option for keeping singletons.
  my @fastq=glob("$$settings{tempdir}/corrected/*.fastq");

  logmsg "Piping to stdout";
  if($isPE){
    @fastq=grep{!/unpaired/} @fastq;
    die "Internal error! More than two read files were found" if(@fastq > 2);
    # Figure out forward/reverse
    my $fwd=(grep{ /_R1_|_1[^\w\d]/} @fastq)[0];
    my $rev=(grep{ /_R2_|_2[^\w\d]/} @fastq)[0];
    system("run_assembly_shuffleReads.pl $fwd $rev");
    die if $?;
  } else {
    system("cat $fastq[0]");
    die if $?;
  }
}

sub usage{
  "$0: Corrects Illumina reads using BayesHammer in SPAdes
  Usage: 
    $0 reads.shuffled.fastq[.gz] > cleaned.shuffled.fastq
    --numcpus 1
    --tempdir <optional>

  "
}
