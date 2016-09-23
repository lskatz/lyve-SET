#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Copy qw/cp/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../../lib";
use LyveSET qw/logmsg @fastaExt @fastqExt/;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help outdir=s numcpus=i)) or die $!;
  
  my($proj)=@ARGV;

  die usage() if(!$proj || $$settings{help});

  # Make the initial directory
  $$settings{outdir}||="$proj/msa/tinsel";
  mkdir $$settings{outdir};
  # Copy over the first two files
  cp("$proj/msa/tree.dnd",$$settings{outdir}) or die "ERROR copying $proj/msa/tree.dnd to $$settings{outdir}: $!";
  cp("$proj/msa/out.pairwiseMatrix.tsv",$$settings{outdir}) or die "ERROR copying $proj/msa/out.pairwiseMatrix.tsv to $$settings{outdir}: $!";
  
  # generate the tips file
  open(METADATA,">","$$settings{outdir}/metadata.tsv") or die "ERROR: could not write to $$settings{outdir}/metadata.tsv: $!";
  print METADATA join("\t",qw(Tip.labels Display.labels))."\n";
  open(FASTA,"$proj/msa/out.aln.fasta") or die "ERROR: could not open $proj/msa/out.aln.fasta: $!";
  while(<FASTA>){
    s/^\s+|\s+$//g;
    next if(!/^>/);
    s/^>//;
    s/\s+.+$//;
    my $sample=$_;

    # remove common suffixes
    $sample=basename($sample,@fastaExt,@fastqExt);

    print METADATA join("\t",$_,$sample)."\n";
  }
  close FASTA;
  close METADATA;

  logmsg "Done! All output files are in $$settings{outdir}";

  return 0;
}

sub usage{
  local $0=basename $0;
  "$0: make the standard files for Tinsel
  Usage: $0 Lyve-SET-project/
  --outdir  proj/msa/tinsel  The output directory
  "
}

