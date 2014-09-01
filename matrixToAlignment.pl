#!/usr/bin/env perl

# Removes uninformative sites from a fasta multiple sequence alignment: Ns, gaps, and invariant sites
# Author: Lee Katz

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;

$0=fileparse $0;
sub logmsg{$|++;print STDERR "$0: @_\n"; $|--;}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help verbose min-distance=i format=s)) or die $!;
  die usage() if($$settings{help});
  $$settings{"min-distance"} ||=0;
  $$settings{format}||="fasta";

  die usage() if($$settings{help});

  my $seqHash=readMatrix($settings);

  return 0;
}

sub readMatrix{
  
}

sub usage{
  "$0: turns a SNP matrix into an alignment
  Usage: $0 matrix.tsv > aln.out.fas
         $0 matrix.tsv | removeUninformativeSites.pl > informative.aln.fas
  --min-distance 0 The distance in bp between allowed SNPs. Does not count sites where all nt are the same.
  --format fasta   The output format
  "
}
