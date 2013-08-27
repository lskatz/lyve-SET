#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use Bio::AlignIO;
use Data::Dumper;
use Getopt::Long;
use threads;
use Thread::Queue;

sub logmsg{$|++;print STDERR "@_\n";$|--;}
exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i));
  die usage() if($$settings{help});
  $$settings{numcpus}||=1;

  my $alignment=$ARGV[0];
  die usage() if(!$alignment);

  logmsg "Reading the aln";
  my %seq;
  my $in=Bio::AlignIO->new(-file=>$alignment);
  while(my $aln=$in->next_aln){
    my @seq=$aln->each_seq;
    for(@seq){
      logmsg $_->id;
      $seq{$_->id}=$_->seq;
    }
  }
  my @seqid=keys(%seq);
  my $numSeq=@seqid;

  # find all distances
  my $printQueue=Thread::Queue->new;
  my $printer=threads->new(\&printer,$printQueue,$settings);
  my $pwQueue=Thread::Queue->new;
  my @thr;
  $thr[$_]=threads->new(\&pairwiseDistanceWorker,\%seq,$pwQueue,$printQueue,$settings) for(0..$$settings{numcpus});
  my @dist;
  # TODO loop based on the seqid
  for(my $i=0;$i<$numSeq;$i++){
    #my $seq1=$seq{$seqid[$i]};
    logmsg "Distances for $seqid[$i]";
    $pwQueue->enqueue($i);
  }

  # close off the threads
  $pwQueue->enqueue(undef) for(@thr);
  $_->join for(@thr);
  $printQueue->enqueue(undef);
  $printer->join;

  return 0;
}
sub pairwiseDistanceWorker{
  my($seqHash,$Q,$printQueue,$settings)=@_;
  my @seqid=keys(%$seqHash);
  my $numSeqs=@seqid;
  while(defined(my $i=$Q->dequeue)){
    my $seq1=$$seqHash{$seqid[$i]};
    for(my $j=$i+1;$j<$numSeqs;$j++){
      my $seq2=$$seqHash{$seqid[$j]};
      my $dist=pairwiseDistance($seq1,$seq2,$settings);
      $printQueue->enqueue(join("\t",$seqid[$i],$seqid[$j],$dist)."\n");
    }
  }
  return 1;
}

sub pairwiseDistance{
  my($seq1,$seq2,$settings)=@_;
  my $length=length($seq1);
  my $pdist=0;
  for(my $i=0;$i<$length;$i++){
    my $nt1=substr($seq1,$i,1);
    my $nt2=substr($seq2,$i,1);
    next if($nt1=~/[N\-\*]/i || $nt2=~/[N\-\*]/i);
    $pdist++ if($nt1 ne $nt2);
  }
  return $pdist;
}
sub printer{
  my($Q,$settings)=@_;
  while(defined(my $thing=$Q->dequeue)){
    print $thing;
  }
  return 1;
}

sub usage{
  "Finds pairwise distances of entries in an alignment file
  Usage: $0 alignment.fasta > pairwise.tsv
  -n numcpus (Default: 1)
  "
}
