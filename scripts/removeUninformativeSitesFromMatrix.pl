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
  my $matrix=readMatrix($settings);
  removeSitesFromMatrix($matrix,$settings);
  #transposeMatrix()
  #print Matrix

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
  return \%matrix_transposed;
}

sub removeSitesFromMatrix{
  my($matrix,$settings)=@_;
  while(my($contig,$posHash)=each(%$matrix)){
    my ($refGenome,$refNt)=each(%$posHash);
    while(my($pos,$genomeHash)=each(%$posHash)){
      # see if it is an invariant site.  If so delete it
      # see if it has a gap
      # see if it has a non-ATCG

    }
  }
}

sub usage{
  "Removes all the uninformative sites in a multiple sequence alignment fasta file: Ns, gaps, and invariant sites
  Usage: $0 < aln.matrix.tsv > informative.matrix.tsv
  -v for verbose (technically makes the script slower)

  Using the following two options will allow you to keep a master MSA list at-hand if you are converting all VCFs in a project
  --gaps-allowed Allow gaps in the alignment
  --ambiguities-allowed Allow ambiguous bases in the alignment
  "
}

