#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/fileparse basename/;
use File::Temp qw/tempdir/;
use File::Copy qw/copy move/;

sub logmsg{local $0=basename($0); print STDERR "@_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i outdir=s regex=s uneven=s));
  $$settings{numcpus}||=1;
  $$settings{outdir}||=basename($0).".out";
  $$settings{regex}||='/^(.*?)(_.*)$/';
  $$settings{uneven}||='die';
  $$settings{uneven}=lc($$settings{uneven});

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
  my $outdir=$$settings{outdir};
  
  # Shuffle all applicable reads
  while(my($basename,$reads)=each(%$pairs)){
    my $outfile="$outdir/$basename.fastq.gz";
    logmsg "Shuffling to $outfile";
    next if(-e $outfile);
    my $tmpfile=shuffle($$reads{1},$$reads{2},$settings);
    system("gzip -v $tmpfile");
    die "ERROR with gzip" if $?;
    system("mv -v $tmpfile.gz $outfile");
    die "ERROR with mv" if $?;
  }
}

sub shuffle{
  my($r1,$r2,$settings)=@_;
  my $tmpdir=tempdir("XXXXXX",TMPDIR=>1,CLEANUP=>1);
  my $shuffled="$tmpdir/shuffled.fastq";
  logmsg "Shuffling to $shuffled";

  open(R1,"<",$r1) or die "ERROR: could not open $r1 for reading: $!";
  open(R2,"<",$r2) or die "ERROR: could not open $r2 for reading: $!";
  open(SHUFFLED,">",$shuffled) or die "ERROR: could not open $shuffled for writing: $!";
  while(<R1>){
    print SHUFFLED $_;           # id1
    $_=<R1>; print SHUFFLED $_;  # seq1
    $_=<R1>; print SHUFFLED $_;  # id1
    $_=<R1>; print SHUFFLED $_;  # qual1

    my $id=<R2>;
    if(!defined($id)){
      my $msg="ERROR: there are not as many reads in $r2 as there are in $r1";
      die $msg if($$settings{uneven} eq 'die');
      logmsg $msg if($$settings{uneven} eq 'warn');
    }
             print SHUFFLED $id; # id2
    $_=<R2>; print SHUFFLED $_;  # seq2
    $_=<R2>; print SHUFFLED $_;  # id2
    $_=<R2>; print SHUFFLED $_;  # qual2
  }

  close R1; close R2; close SHUFFLED;
  return $shuffled;
}

sub usage{
  local $0=fileparse $0;
  "$0: shuffle a set of fastq files into a directory
  Usage: $0 *.fastq.gz -o outdir
  --outdir  output directory
  --regex   A regular expression with parentheses to capture the basename of both F/R files. Also must have parentheses to capture the second half of the filename.
  --uneven  How to deal with uneven reads.  Options: die, warn
  "
}
