#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/fileparse basename/;
use File::Copy qw/copy move/;
use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i outdir=s regex=s)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{outdir}||=basename($0).".out";
  $$settings{regex}||='/^(.*?)(_.*)$/';

  my @splitFastq=@ARGV;
  die "ERROR: no fastq files were given!\n".usage() if(!@splitFastq);
  
  mkdir($$settings{outdir}) if(!-d $$settings{outdir});

  my $pairs=readFastqs(\@splitFastq,$settings);
  shufflePairs($pairs,$settings);
}

sub readFastqs{
  my($fastq,$settings)=@_;
  
  # Assign pairs based on the fastq regex
  my $regex=$$settings{regex};
  $regex=~s/^\/|\/$//g; # remove slashes in the regex
  my %fastqPair;
  for my $filename(@$fastq){
    $filename=~/$regex/;
    my($basename,$therest)=($1,$2);
    my $readNumber;
    if($therest=~/_R([12])_/){
      $readNumber=$1;
    } elsif($therest=~/_([12])\.f/){
      $readNumber=$1;
    } else {
      die "ERROR: could not parse the read number from file $filename in the second part of the file =>$therest<=";
    }

    die "ERROR: trying to set $filename as read number $readNumber for $basename, but it already exists as ".$fastqPair{$basename}{$readNumber} if($fastqPair{$basename}{$readNumber});
    $fastqPair{$basename}{$readNumber}=$filename;
  }

  return \%fastqPair;
}

sub shufflePairs{
  my($pairs,$settings)=@_;

  my $Q=Thread::Queue->new;

  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&shuffleWorker,$Q,$settings);
  }

  # Shuffle all applicable reads
  while(my($basename,$reads)=each(%$pairs)){
    $Q->enqueue([$basename,$reads]);
  }

  # Join the threads
  $Q->enqueue(undef) for(@thr);
  $_->join for(@thr);
}

sub shuffleWorker{
  my($Q,$settings)=@_;
  my $outdir=$$settings{outdir};

  while(defined(my $tmp=$Q->dequeue)){
    my($basename,$reads)=@$tmp;

    logmsg "Looking for $basename";
    my $outfile="$outdir/$basename.fastq.gz";
    next if(-e $outfile);

    system("run_assembly_shuffleReads.pl $$reads{1} $$reads{2} | gzip -c > $outfile.tmp");
    die "ERROR with run_assembly_shuffleReads.pl" if $?;
    system("mv -v $outfile.tmp $outfile");
  }
}

sub usage{
  local $0=fileparse $0;
  "$0: shuffle a set of fastq files into a directory. Uses run_assembly_shuffleReads.pl.
  Usage: $0 *.fastq.gz -o outdir
  --outdir  output directory
  --regex   A regular expression with parentheses to capture the basename of both F/R files. Also must have parentheses to capture the second half of the filename.
  "
}
