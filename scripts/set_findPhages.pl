#!/usr/bin/env perl
# Author:Lee Katz <lkatz@cdc.gov>
# Thanks: Darlene Wagner for giving me this idea

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use File::Temp qw/tempdir/;
use List::Util qw/min max/;
use List::MoreUtils qw(uniq);
use Bio::Perl;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;
use lib "$FindBin::RealBin/../lib/lib/perl5";
use Number::Range;

local $0=fileparse $0;
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i tempdir=s));
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("phastXXXXXX",CLEANUP=>1,TMPDIR=>1);

  my $fasta=$ARGV[0];
  die usage() if(!$fasta || $$settings{help});

  logmsg "Tempdir is $$settings{tempdir}";
  my $ranges=phast($fasta,$settings);

  # Print the ranges to stdout
  for my $r(@$ranges){
    print join("\t",@$r)."\n";
  }
  return 0;
}

sub phast{
  my($fasta,$settings)=@_;

  my $tempdir=tempdir("$$settings{tempdir}/phastXXXXXX",CLEANUP=>1,PREFIX=>$$settings{tempdir});
  my $db="$FindBin::RealBin/../lib/phast/phast.faa";
  logmsg "Running blastx against $db";

  # Better parallelization: one fasta entry per cpu.
  # Split the query into multiple files and then figure out
  # how many cpus per blast job we need.
  my $i=0;
  my $seqin=Bio::SeqIO->new(-file=>$fasta);
  while(my $seq=$seqin->next_seq){
    my $file="$tempdir/$i.fna";
    open(SEQOUT,">",$file) or die "ERROR: could not write seq to temp file $file: $!";
    print SEQOUT join("\n",">".$seq->id,$seq->seq,"");
    close SEQOUT;
    $i++;
  }
  my $threadsPerBlast=int($$settings{numcpus}/$i);
  $threadsPerBlast=1 if($threadsPerBlast<1);

  # Perform blast on these split files.
  system("ls $tempdir/*.fna | xargs -I {} -P $$settings{numcpus} -n 1 blastx -query {} -db $db -evalue 0.05 -outfmt 6 -num_threads $threadsPerBlast -out {}.bls");
  die "ERROR with blastx: $!" if $?;
  #my $allResults=`blastx -query '$fasta' -db $db -evalue 0.05 -outfmt 6 -num_threads $$settings{numcpus}`;
  my $allResults=`cat $tempdir/*.fna.bls`;
  die "ERROR with cat on $tempdir/*.fna.bls" if($? || !$allResults);

  logmsg "Parsing results";
  my(%range);
  for my $result(split(/\n/,$allResults)){
    $result=~s/^\s+|\s+$//g; # trim
    my ($contig,$hit,$identity,$length,$gaps,$mismatches,$sstart,$send,$qstart,$qend,$e,$score)=split /\t/, $result;
    next if($score < 50 || $length < 20);
    
    # Add these coordinates to ranges
    $range{$contig}||=Number::Range->new;
    my $lo=min($sstart,$send);
    my $hi=max($sstart,$send);
    
    no warnings;
    $range{$contig}->addrange($lo..$hi);
  }

  # Translate the ranges found in the Range objects into 
  # an array of [contig,start,stop]
  my @range;
  while(my($contig,$rangeObj)=each(%range)){
    my $rangeStr=$rangeObj->range;
    while($rangeStr=~/(\d+)\.\.(\d+),?/g){
      push(@range,[$contig,$1,$2]);
    }
  }

  return \@range;
}

sub readFasta{
  my($fasta,$settings)=@_;
  my $in=Bio::SeqIO->new(-file=>$fasta);
  my %seq;
  while(my $seq=$in->next_seq){
    $seq{$seq->id}=$seq->seq;
  }
  return %seq;
}

sub usage{
  "Finds phages in a fasta file using phast and PhiSpy
  Usage: $0 file.fasta
  --numcpus 1
  --tempdir tmp/
  "
}
