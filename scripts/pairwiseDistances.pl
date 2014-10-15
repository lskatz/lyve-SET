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
  GetOptions($settings,qw(help numcpus=i estimate samplingFrequency=s));
  die usage() if($$settings{help});
  $$settings{numcpus}||=1;
  $$settings{samplingFrequency}||=0.25;

  # TODO change this so that you can read stdin
  my $alignment=$ARGV[0];
  die usage() if(!$alignment);

  my %seq;
  my $in=Bio::AlignIO->new(-file=>$alignment,-idlength=>30);
  while(my $aln=$in->next_aln){
    my @seq=$aln->each_seq;
    for(@seq){
      $seq{$_->id}=lc($_->seq);
    }
  }
  my @seqid=keys(%seq);
  my $numSeq=@seqid;
  logmsg "Loaded $numSeq sequences for comparison";

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

  # get the status while the threads are still working on it
  while($pwQueue->pending > 0){
    for(0..28){ # leave the last second for the end just in case the queue gets added onto
      sleep 1;
      last if($pwQueue->pending == 0);
      logmsg $pwQueue->pending." queries left in the queue" if($_ == 0);
    }
    sleep 1; # just in case something gets added onto the print queue
  }
  logmsg "No more left in the queue";

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
      my($s1,$s2)=sort {$a cmp $b} ($seqid[$i],$seqid[$j]); # to have the pairwise in some order
      $printQueue->enqueue(join("\t",$s1,$s2,$dist)."\n");
    }
  }
  return 1;
}

sub pairwiseDistance{
  my($seq1,$seq2,$settings)=@_;
  #$seq1=lc($seq1); $seq2=lc($seq2); # I added this when reading the sequences
  my $length=length($seq1);
  my $pdist=0;
  for(my $i=0;$i<$length;$i++){
    next if($$settings{estimate} && rand() >$$settings{samplingFrequency});
    my $nt1=substr($seq1,$i,1);
    my $nt2=substr($seq2,$i,1);
    next if($nt1=~/[^atcg]/ || $nt2=~/[^atcg]/ || $nt1 eq $nt2);
    $pdist++;
  }
  $pdist=int($pdist / $$settings{samplingFrequency}) if($$settings{estimate});
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
  -e Estimate the number of pairwise distances using random sampling. 1/4 of all pairwise bases will be analyzed instead of 100%.
  -s 0.25 (to be used with -e) The frequency at which to analyze positions for pairwise differences
  "
}
