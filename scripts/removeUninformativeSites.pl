#!/usr/bin/env perl

# Removes uninformative sites from a fasta multiple sequence alignment: Ns, gaps, and invariant sites
# Author: Lee Katz

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO;

sub logmsg{$|++;print STDERR "@_\n"; $|--;}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help verbose ambiguities-allowed gaps-allowed sort=s)) or die $!;
  die usage() if($$settings{help});
  $$settings{"ambiguities-allowed"} ||=0;
  $$settings{"gaps-allowed"} ||=0;
  $$settings{sort}||="";
  $$settings{sort}=lc($$settings{sort});

  ## read in the fasta file into @seq and %seq, and keep the deflines
  my $in;
  if(defined($ARGV[0]) && -f $ARGV[0]){
    $in=Bio::SeqIO->new(-file=>$ARGV[0]);
  } else {
    $in=Bio::SeqIO->new(-fh=>\*STDIN,-format=>"fasta");
  }
  my($length,$defline,@defline,@seq);
  while(my $seqObj=$in->next_seq){
    push(@seq,$seqObj->seq);
    push(@defline,$seqObj->id);
    $length=$seqObj->length;
  }
  $in->close;

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
  my %seq;
  @seq{@defline}=@seq;
  # sort the sequences by identifier
  my @sortedId;
  if($$settings{sort} eq 'alpha'){
    @sortedId=sort {$a cmp $b} @defline;
  } elsif($$settings{sort} eq 'num'){
    @sortedId=sort {$a <=> $b} @defline;
  } else {
    @sortedId=@defline;
  }
  for (my $i=0;$i<@sortedId;$i++){
    my $id=$sortedId[$i];
    my $sequence=$seq{$id};
    print ">$id\n";
    for my $pos(@pos){
      print substr($sequence,$pos,1);
    }
    print "\n";
  }
  return 0;


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
  --sort alpha,num Sort the sequences by their deflines.  Values for sort are either ALPHA or NUM (case insensitive)
  "
}

