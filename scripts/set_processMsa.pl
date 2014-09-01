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

sub logmsg {local $0=basename $0;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}
exit main();

sub main{
  local $0=basename $0;
  my $settings={};
  GetOptions($settings,qw(auto groups=s@ treePrefix=s alnPrefix=s pairwisePrefix=s fstPrefix=s eigenPrefix=s numcpus=i msaDir=s tempdir=s force help)) or die $!;
  $$settings{treePrefix}||="";
  $$settings{alnPrefix}||="";
  $$settings{pairwisePrefix}||="";
  $$settings{fstPrefix}||="";
  $$settings{eigenPrefix}||="";
  $$settings{numcpus}||=1;
  $$settings{auto}||=0;
  $$settings{msaDir}||=0;
  $$settings{force}||=0;
  die usage() if($$settings{help});

  if($$settings{auto} || $$settings{msaDir}){
    if($$settings{auto}){
      $$settings{msaDir}="msa" if(-e "msa/out.aln.fas");
      $$settings{msaDir}="." if(-e "./out.aln.fas");
      die "ERROR: --auto was set but I could not find out.aln.fas in either this directory or in ./msa/\n".usage() if(!$$settings{msaDir});
    }
    $$settings{auto}=1;  # explicitly set auto if auto or msaDir is set

    my $dir="$$settings{msaDir}"; # save me on typing the next few lines
    $$settings{treePrefix}||="$dir/RAxML";
    $$settings{alnPrefix}||="$dir/informative";
    $$settings{pairwisePrefix}||="$dir/pairwise";
    $$settings{fstPrefix}||="$dir/fst";
    $$settings{eigenPrefix}||="$dir/eigen";
  }
  $$settings{tempdir}||="$$settings{msaDir}/tmp";
  mkdir $$settings{tempdir} if(!-d $$settings{tempdir});

  my $infile;
  if(@ARGV){
    $infile=$ARGV[0];
  } elsif($$settings{auto}){
    $infile="$$settings{msaDir}/out.aln.fas";
  } else {
    die "ERROR: need input alignment file\n".usage();
  }
  if(!-f $infile){
    die "ERROR: could not find file $infile\n".usage();
  }
  logmsg "Input alignment file is ".File::Spec->rel2abs($infile);

  # Start off the threads
  my @thr;
  $thr[0]=threads->new(\&distanceStuff,$infile,$settings);
  $thr[1]=threads->new(\&phylogenies,$infile,$settings);
  $_->join for(@thr); @thr=();

  # Things that depend on a tree and pairwise output
  Fst($infile,$$settings{pairwisePrefix},$$settings{treePrefix},$$settings{fstPrefix},$settings) if($$settings{fstPrefix} && $$settings{treePrefix} && $$settings{pairwisePrefix});

  #rmdir $$settings{tempdir};  # don't force this rmdir in case it contains files. This script should remove all tmp files before exiting.

  return 0;
}

####
## distance metrics stuff
####

sub distanceStuff{
  my($infile,$settings)=@_;
  logmsg "Calculating distances";

  my $pairwise=pairwiseDistance($infile,$$settings{pairwisePrefix},$settings) if($$settings{pairwisePrefix});
  eigen($pairwise,$$settings{eigenPrefix},$settings) if($$settings{eigenPrefix} && $pairwise);
}

sub pairwiseDistance{
  my($infile,$prefix,$settings)=@_;
  my $outfile="$prefix.tsv";
  my $matrix="$prefix.matrix.tsv";
  if(-f $outfile && !$$settings{force}){
    logmsg "$outfile was found. I will not perform pairwise distances again without --force";
    return $outfile;
  }
  logmsg "Calculating pairwise distances";
  system("pairwiseDistances.pl --numcpus $$settings{numcpus} '$infile' | sort -k3,3n > '$outfile'");
  die if $?;
  system("pairwiseTo2d.pl < '$outfile' > '$matrix'");
  die if $?;
  
  return $outfile;
  # TODO inter and intra group distances
}

sub eigen{
  my($pairwise,$prefix,$settings)=@_;
  if(-f "$prefix.tsv" && !$$settings{force}){
    logmsg "The eigen vector file was found in $prefix.tsv. I will not recreate it without --force";
    return "$prefix.tsv";
  }
  system("set_indexCase.pl $pairwise | sort -k2,2nr > $prefix.tsv");
  logmsg "ERROR in set_indexCase.pl: $!" if $?;
  return "$prefix.tsv";
}


######
## phylogeny subroutines
#####


sub phylogenies{
  my($inAln,$settings)=@_;
  logmsg "Calculating phylogenies";

  my($informativeAln,$tree);
  $informativeAln=removeUninformativeSites($inAln,$$settings{alnPrefix},$settings) if($$settings{alnPrefix});
  $tree=inferPhylogeny($informativeAln,$$settings{treePrefix},$settings) if($$settings{treePrefix});
  return($informativeAln,$tree);
}

sub removeUninformativeSites{
  my($inAln,$outPrefix,$settings)=@_;
  my $informative="$outPrefix.aln.fas";
  if(-f $informative && !$$settings{force}){
    logmsg "$informative was found.  I will not recalculate without --force.";
    return $informative;
  }
  logmsg "Removing uninformative sites from the alignment and putting it into $informative";
  logmsg "  removeUninformativeSites.pl < '$inAln' > '$informative'";
  system("removeUninformativeSites.pl < '$inAln' > '$informative'");
  die if $?;
  return $informative;
}

# TODO put back in phyml
sub inferPhylogeny{
  my($inAln,$prefix,$settings)=@_;
  my $treeFile="$prefix.RAxML_bipartitions";
  if(-f "$prefix.RAxML_bipartitions" && !$$settings{force}){
    logmsg "$prefix.RAxML_bipartitions was found. I will not recalculate the phylogeny without --force.";
    return $treeFile;
  }
  
  my $numTaxa=`grep -c '>' '$inAln'`; die if $?;
  chomp($numTaxa);
  if($numTaxa < 4){
    logmsg "Only $numTaxa in the alignment. I will not determine the phylogeny";
    return "";
  }

  system("rm -fv $$settings{tempdir}/RAxML*");
  logmsg "Running raxml";
  my $alnInAbs=File::Spec->rel2abs($inAln);
  my $command="cd $$settings{tempdir}; launch_raxml.sh $alnInAbs suffix";
  logmsg "  $command";
  system($command);
  die if $?;

  # Move those files over when finished
  for (qw(RAxML_bestTree RAxML_bipartitionsBranchLabels RAxML_bipartitions RAxML_bootstrap RAxML_info)){
    system("mv -v $$settings{tempdir}/$_.suffix $prefix.$_");
    die "ERROR: could not move $$settings{tempdir}/$_.suffix to $prefix.$_: $!" if $?;
  }
  return $treeFile
}

#### 
# things that depend on pairwise and trees
####

sub Fst{
  my($inAln,$pairwisePrefix,$treePrefix,$fstPrefix,$settings)=@_;
  my $fstTree="$fstPrefix.fst.dnd";
  if(-f $fstTree && !$$settings{force}){
    logmsg "$fstTree fst tree was found.  Not recalculating without --force.";
    return $fstTree;
  }
  if(!-f "$treePrefix.RAxML_bipartitions"){
    logmsg "Tree was not created. Will not perform Fst on an empty tree";
    return "";
  }
  # TODO multithread this?
  system("applyFstToTree.pl --numcpus $$settings{numcpus} -t $treePrefix.RAxML_bipartitions -p $pairwisePrefix.tsv --outprefix $fstPrefix --outputType averages > $fstPrefix.avg.tsv");
  die if $?;
  system("applyFstToTree.pl --numcpus $$settings{numcpus} -t $treePrefix.RAxML_bipartitions -p $pairwisePrefix.tsv --outprefix $fstPrefix --outputType samples > $fstPrefix.samples.tsv");
  die if $?;
  return $fstTree;
}

sub usage{
  local $0=basename $0;
  "Process an MSA from a hqSNP pipeline and get useful information
  Usage: $0 file.fasta
  #-g   groups.txt       Text files of genome names, one per line (multiple -g are encouraged).
                        Groups will help with Fst and with pairwise distances.
  -n   numcpus
  --force               Files will be overwritten even if they exist
  --tempdir tmp/        Place to put temporary files

OUTPUT
  -tre treePrefix
  -aln informative.aln  Informative alignment, created by removeUninformativeSites.pl
  -p   pairwisePrefix   Pairwise distances
  -fst fstPrefix        Fixation index output files
  -e   eigenPrefix      Eigenvalue output files (connectedness)

  --auto                Indicate that you are currently in a SET directory and that you would like all outputs, overwriting any other output parameters
  --msaDir              Indicate a directory to perform everything in. --auto will be set if --msaDir is invoked
  "
}
