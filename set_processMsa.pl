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

my $sge=Schedule::SGELK;
sub logmsg {local $0=basename $0;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}
exit main();

sub main{
  local $0=basename $0;
  my $settings={};
  GetOptions($settings,qw(auto groups=s@ treePrefix=s alnPrefix=s pairwisePrefix=s fstFile=s eigenPrefix=s numcpus=i)) or die $!;
  my $infile=$ARGV[0];
  die "ERROR: need input alignment file\n".usage() if(!$infile || !-f $infile);
  $$settings{treePrefix}||="";
  $$settings{alnPrefix}||="";
  $$settings{pairwisePrefix}||="";
  $$settings{fstFile}||="";
  $$settings{eigenPrefix}||="";
  $$settings{numcpus}||=1;

  if($$settings{auto}){
    $$settings{treePrefix}||="msa/RAxML.bipartitions.informative";
    $$settings{alnPrefix}||="msa/informative";
    $$settings{pairwisePrefix}||="msa/pairwise";
    $$settings{fstFile}||="msa/fst.tsv";
    $$settings{eigenPrefix}||="msa/eigen";
  }

  # Kick off pairwise distances, Fst, eigen stuff
  distanceStuff($infile,$settings);
  # Kick off tree stuff
  phylogenies($infile,$settings);
  # Things that depend on a tree and pairwise output
  $sge->wrapItUp(); # finish all prereqs
  Fst($$settings{fstFile},$settings) if($$settings{fstFile});

  return 0;
}

sub distanceStuff{
  my($infile,$settings)=@_;

  my($pairwise,$fst,$eigen);

  $pairwise=pairwiseDistance($infile,$$settings{pairwisePrefix},$settings) if($$settings{pairwisePrefix});

  # Eigen depends on pairwise, so wait for it to finish.
  $sge->wrapItUp();
  eigen($pairwise,$$settinsg{eigenPrefix},$settings) if($$settings{eigenPrefix} && $pairwise);
}

sub pairwiseDistance{
  my($infile,$prefix,$settings)=@_;
  my $outfile="$prefix.tsv";
  $sge->pleaseExecute("pairwiseDistances.pl --numcpus $$settings{numcpus} '$infile' > '$outfile'",{numcpus=>$$settings{numcpus},jobname=>"pairwiseDistance");
  die if $?;
  
  return $outfile;
  # TODO inter and intra group distances
}

sub eigen{
  my($infile,$prefix,$settings)=@_;
  logmsg "TODO make eigenvector script";
}


sub phylogenies{
  my($inAln,$settings)=@_;

  my $informativeAln=$inAln; # if an informative MSA is not specified, then this one will do.
  $informativeAln=removeUninformativeSites($inAln,$$settings{alnPrefix},$settings) if($$settings{alnPrefix});
  $tree=inferPhylogeny($informativeAln,$$settings{treePrefix},$settings) if($$settings{treePrefix});
}

sub removeUninformativeSites{
  my($inAln,$outPrefix,$settings)=@_;
  $sge->pleaseExecute("removeUninformativeSites.pl < '$inAln' > '$informative'",{jobname=>"removeUninformativeSites",numcpus=>1});
  die if $?;
  return $informative;
}

sub phylogenies{
  my($inAln,$prefix,$settings)=@_;
  
}

sub Fst{
  my($inAln,$inTree,$$settings{fstFile},$settings)=@_;
}

sub usage{
  local $0=basename $0;
  "Process an MSA from a hqSNP pipeline and get useful information
  Usage: $0 file.fasta
  -g   groups.txt       Text files of genome names, one per line (multiple -g are encouraged).
                        Groups will help with Fst and with pairwise distances.
  -n   numcpus
  Output options:
  -t   treePrefix
  -aln informative.aln  Informative alignment, created by removeUninformativeSites.pl
  -p   pairwisePrefix   Pairwise distances
  -f   fstPrefix        Fixation index output files
  -e   eigenPrefix      Eigenvalue output files
  --auto                Indicate that you are currently in a SET directory and that you would like all outputs, overwriting any other output parameters
  "
}
