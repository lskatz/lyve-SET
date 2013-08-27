#!/usr/bin/env perl

# Removes uninformative sites from a fasta multiple sequence alignment: Ns, gaps, and invariant sites
# Author: Lee Katz

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

sub logmsg{$|++;print STDERR "@_\n"; $|--;}
my $settings={};
GetOptions($settings,qw(help verbose));
die usage() if($$settings{help});

## read in the fasta file into @seq and %seq, and keep the deflines
my($length,$defline,@defline,@seq,%seq);
while(<>){
  chomp; 
  if(/^>/){
    s/>//; 
    $defline=$_;
    push(@defline,$defline);
    next;
  } 
  $seq{$defline}=$_; 
  push(@seq,$_); 
  $length=length($_);
} 

## read informative positions into @pos
# compare each sequence to the reference sequence (the first sequence)
my $refSeq=shift(@seq); 
my $refId=shift(@defline); 
my (%aln,@pos);
my $numOtherSeq=@seq;
my $informativeCount=0;
POSITION:for(my $j=0;$j<$length;$j++){ 
  my $informative=0; 
  for(my $i=0;$i<$numOtherSeq;$i++){
    my $sequence=$seq[$i];
    my $nt=substr($sequence,$j,1); 
    $informative=1 if($nt ne substr($refSeq,$j,1)); 
    next POSITION if($nt=~/[nN\-]/); 
  } 
  next if(!$informative);
  my $refNt=substr($refSeq,$j,1);
  next if($refNt=~/[nN\-]/); # if it's informative, but the ref base isn't good, then skip
  push(@pos,$j);
  
  logmsg $j if($informativeCount % 100 == 0 && $$settings{verbose});
  $informativeCount++;
}

## print all nucleotides found at informative positions
# bring back the reference sequence
unshift(@seq,$refSeq);
unshift(@defline,$refId);
for (my $i=0;$i<@seq;$i++){
  my $sequence=$seq[$i];
  my $id=$defline[$i];
  print ">$id\n";
  for my $pos(@pos){
    print substr($sequence,$pos,1);
  }
  print "\n";
}

return 0;

sub usage{
  "Removes all the uninformative sites in a multiple sequence alignment fasta file: Ns, gaps, and invariant sites
  Usage: $0 < aln.fasta > informative.fasta
  -v for verbose (technically makes the script slower)
  "
}
