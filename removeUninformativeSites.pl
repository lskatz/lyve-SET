#!/usr/bin/env perl

# Removes uninformative sites from a fasta multiple sequence alignment: Ns, gaps, and invariant sites
# Author: Lee Katz

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

sub logmsg{$|++;print STDERR "@_\n"; $|--;}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help verbose ambiguities-allowed gaps-allowed));
  die usage() if($$settings{help});
  $$settings{"ambiguities-allowed"} ||=0;
  $$settings{"gaps-allowed"} ||=0;

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
  my $removeAmbiguities=!$$settings{"ambiguities-allowed"};
  my $removeGaps=!$$settings{"gaps-allowed"};
  my $refSeq=shift(@seq); 
  my $refId=shift(@defline); 
  my (%aln,@pos);
  my $numOtherSeq=@seq;
  my $informativeCount=0;
  POSITION:for(my $j=0;$j<$length;$j++){ 
    my $refNt=substr($refSeq,$j,1);
    next if($removeAmbiguities && uc($refNt) eq 'N'); # if it's informative, but the ref base isn't good, then skip
    next if($removeGaps && $refNt eq '-');            # if it's informative, but the ref base isn't good, then skip

    # see if the rest of this column shows that it is informative
    my $informative=0; # guilty until proven innocent
    for(my $i=0;$i<$numOtherSeq;$i++){
      #my $sequence=$seq[$i];
      my $nt=substr($seq[$i],$j,1); 
      die "ERROR: Sequence $i does not have a nucleotide at position $j! Is the MSA flush?" if(!$nt);
      $informative=1 if($nt ne $refNt);  # It's informative, but you have to continue reviewing
                                         # the other sequences if you are ignoring N or gap columns.
      next POSITION if($removeAmbiguities && uc($nt) eq 'N'); # no chance of redemption here
      next POSITION if($removeGaps && $nt eq '-');            # no chance of redemption here
    } 
    next if(!$informative);
    push(@pos,$j); # retain this informative position
    
    $informativeCount++;
    logmsg $j if($$settings{verbose});
  }

  ## Print all nucleotides found at informative positions.
  # Bring back the reference sequence.
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
}

sub usage{
  "Removes all the uninformative sites in a multiple sequence alignment fasta file: Ns, gaps, and invariant sites
  Usage: $0 < aln.fasta > informative.fasta
  -v for verbose (technically makes the script slower)

  Using the following two options will allow you to keep a master MSA list at-hand if you are converting all VCFs in a project
  --gaps-allowed Allow gaps in the alignment
  --ambiguities-allowed Allow ambiguous bases in the alignment
  "
}

