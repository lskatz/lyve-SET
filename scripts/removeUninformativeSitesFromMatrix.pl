#!/usr/bin/env perl

# Removes uninformative sites from a matrix: Ns, gaps, and invariant sites
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

  ## read in the matrix into a hash of hashes
  my ($posArr,$matrix)=readMatrix($settings);
  removeSitesFromMatrix($matrix,$settings);

  # first row needs to be the loci; next rows are genomes/nts
  my $matrix_transposed=transposeMatrix($matrix,$settings);
  printMatrix($matrix_transposed,$posArr,$settings);

  return 0;
}

sub readMatrix{
  my($settings)=@_;
  my $MATRIX; # file pointer
  my %matrix; # hash of nts/positions
  my %matrix_transposed;
  if(defined($ARGV[0]) && -e $ARGV[0]){
    open($MATRIX,$ARGV[0]) or die "ERROR: could not open matrix $ARGV[0]: $!";
  } else {
    $MATRIX=\*STDIN;
  }
  my $header=<$MATRIX>; chomp $header;
  my(undef,@pos)=split(/\t/,$header);
  my $numPos=@pos;
  while(<$MATRIX>){
    chomp;
    my ($genome,@nt)=split(/\t/,$_);
    for (my $i=0;$i<$numPos;$i++){
      my($contig,$pos)=split(/:/,$pos[$i]);
      $matrix{$genome}{$contig}{$pos}=$nt[$i];
      $matrix_transposed{$contig}{$pos}{$genome}=$nt[$i];
    }
  }

  return (\@pos,\%matrix_transposed);
}

sub removeSitesFromMatrix{
  my($matrix,$settings)=@_;
  while(my($contig,$posHash)=each(%$matrix)){
    my ($refGenome,$refNt)=each(%$posHash);
    my $REFNT=uc($refNt);
    while(my($pos,$genomeHash)=each(%$posHash)){
      my ($is_variant,$is_indel,$is_ambiguous)=(0,0,0);
      for my $nt(values(%$genomeHash)){
        my $NT = uc($nt);
        # see if it is an invariant site.
        $is_variant=1 if($NT ne $REFNT);
        # see if it has a gap
        $is_indel=1 if($NT eq '*' || $NT eq '-' || length($NT) > 1 || $REFNT eq '*' || $REFNT eq '-' || length($REFNT) > 1);
        # see if it has a non-ATCG
        $is_ambiguous=1 if($NT=~/^\w+$/ && $NT !~/^[ATGC]$/);
      }

      # Remove or keep the NT based on the boolean logic and the settings requested
      delete($$matrix{$contig}{$pos}) if(!$$settings{'gaps-allowed'} && $is_indel);
      delete($$matrix{$contig}{$pos}) if(!$$settings{'ambiguities-allowed'} && $is_ambiguous);
      delete($$matrix{$contig}{$pos}) if(!$is_variant);

      if($$settings{verbose}){
        logmsg "Deleted $contig:$pos" if(!$$matrix{$contig}{$pos});
        logmsg "Kept $contig:$pos" if($$matrix{$contig}{$pos});
      }
    }
  }
}

sub transposeMatrix{
  my($matrix,$settings)=@_;
  my %transposed;
  while(my($contig,$posHash)=each(%$matrix)){
    while(my($pos,$genomeHash)=each(%$posHash)){
      while(my($genome,$nt)=each(%$genomeHash)){
        $transposed{$genome}{$contig}{$pos}=$nt;
      }
    }
  }
  return \%transposed;
}

sub printMatrix{
  my($matrix,$posArr,$settings)=@_;
  
  # get a list of contigs and hashes
  my %contig;
  for my $posKey(@$posArr){
    my($contig,$pos)=split(/:/,$posKey);
    push(@{ $contig{$contig} },$pos);
  }
  #sort positions
  for (values(%contig)){
    my @sorted=sort {$a<=>$b} @$_;
    $_=\@sorted;
  }

  # Sort the header positional keys
  my @sortedPosKey=sort {
    my ($contigA,$posA)=split(/:/,$a);
    my ($contigB,$posB)=split(/:/,$b);
    return $contigA cmp $contigB if($contigA ne $contigB);
    return $posA <=> $posB;
  } @$posArr;
  print join("\t",".",@sortedPosKey)."\n";

  # Body of matrix:
  # Do a sorted print.
  while(my($genome,$contigHash)=each(%$matrix)){
    for my $contig(sort {$a cmp $b} keys(%$contigHash)){
      my $posHash=$$matrix{$genome}{$contig};
      my $line="$genome\t";
      for my $pos(sort {$a <=> $b} keys(%$posHash)){
        $line.=$$posHash{$pos}."\t";
      }
      $line=~s/\t+$//;
      print $line."\n";
    }
  }

  return 1;
}

sub usage{
  "Removes all the uninformative sites in a multiple sequence alignment fasta file: Ns, gaps, and invariant sites
  Usage: $0 < aln.matrix.tsv > informative.matrix.tsv
  -v for verbose

  --gaps-allowed Allow gaps in the matrix
  --ambiguities-allowed Allow ambiguous bases in the matrix (non-ATCG)
  "
}

