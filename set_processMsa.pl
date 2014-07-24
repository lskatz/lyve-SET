#!/usr/bin/env perl
# Process a resulting MSA: remove uninformative sites; find pairwise distances;
# find Fst; make a tree; calculate the eigenvector
# Author: Lee Katz <lkatz@cdc.gov>

use FindBin;
use lib "$FindBin::RealBin/lib";

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Spec;
use threads;
use Thread::Queue;
use Schedule::SGELK;

sub logmsg {local $0=basename $0;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}
exit main();

sub main{
  local $0=basename $0;
  my $settings={};
  GetOptions($settings,qw(output=s treePrefix=s alignmentInformative=s pairwiseFile=s fstFile=s eigenPrefix=s)) or die $!;
  my $infile=$ARGV[0];
  die "ERROR: need input alignment file\n".usage() if(!$infile || !-f $infile);
  $$settings{treePrefix}||="$0.tree";
  $$settings{alignmentInformative}||="$0.informative";
  $$settings{pairwiseFile}||="$0.pairwise";
  $$settings{fstFile}||="$0.fst";
  $$settings{eigenPrefix}||="$0.eigen";


  # Kick off pairwise distances, Fst, eigen stuff
  distanceStuff($infile,$settings);
  # Kick off tree stuff
  phylogenies($infile,$settings);
}

sub usage{
  local $0=basename $0;
  "Process an MSA from a hqSNP pipeline and get useful information
  Usage: $0 file.fasta
  Output options:
  -t treePrefix
  -o output.aln       Output alignment
  -a informative.aln  Informative alignment, created by removeUninformativeSites.pl
  -p pairwise.tsv     Pairwise distances
  -f fstPrefix        Fixation index output files
  -e eigenPrefix      Eigenvalue output files
  "
}
