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
  $$settings{tempdir}||="tmp";
  mkdir $$settings{tempdir} if(!-d $$settings{tempdir});
  die usage() if($$settings{help});

  if($$settings{auto} || $$settings{msaDir}){
    if($$settings{auto}){
      $$settings{msaDir}="msa" if(-e "msa/out.aln.fas");
      $$settings{msaDir}="." if(-e "./out.aln.fas");
      die "ERROR: --auto was set but I could not find out.aln.fas in either this directory or in ./msa/\n".usage() if(!$$settings{msaDir});
    }
    my $dir="$$settings{msaDir}"; # save me on typing
    $$settings{treePrefix}||="$dir/RAxML";
    $$settings{alnPrefix}||="$dir/informative";
    $$settings{pairwisePrefix}||="$dir/pairwise";
    $$settings{fstPrefix}||="$dir/fst";
    $$settings{eigenPrefix}||="$dir/eigen";

    $$settings{auto}=1;  # explicitly set auto
  }

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

  logmsg "Calculating distances";
  distanceStuff($infile,$settings);
  logmsg "Calculating phylogenies";
  phylogenies($infile,$settings);
  # Things that depend on a tree and pairwise output
  Fst($infile,$$settings{pairwisePrefix},$$settings{treePrefix},$$settings{fstPrefix},$settings) if($$settings{fstPrefix} && $$settings{treePrefix} && $$settings{pairwisePrefix});

  return 0;
}

####
## distance metrics stuff
####

sub distanceStuff{
  my($infile,$settings)=@_;

  my($pairwise,$fst,$eigen);
  my $sge=Schedule::SGELK->new(keep=>1,numcpus=>$$settings{numcpus});

  $pairwise=pairwiseDistance($infile,$$settings{pairwisePrefix},$sge,$settings) if($$settings{pairwisePrefix});

  # Eigen depends on pairwise, so wait for it to finish.
  $sge->wrapItUp();
  eigen($pairwise,$$settings{eigenPrefix},$sge,$settings) if($$settings{eigenPrefix} && $pairwise);

  $sge->wrapItUp();  # final wrap-up before leaving this sub
}

sub pairwiseDistance{
  my($infile,$prefix,$sge,$settings)=@_;
  my $outfile="$prefix.tsv";
  if(-f $outfile){
    logmsg "$outfile was found. I will not perform pairwise distances again";
    return $outfile;
  }
  logmsg "Calculating pairwise distances";
  $sge->pleaseExecute("pairwiseDistances.pl --numcpus $$settings{numcpus} '$infile' | sort -k3,3n > '$outfile'",{numcpus=>$$settings{numcpus},jobname=>"pairwiseDistance"});
  
  return $outfile;
  # TODO inter and intra group distances
}

# TODO put this in a separate script
sub eigen{
  my($pairwise,$prefix,$sge,$settings)=@_;
  $sge->pleaseExecute("set_indexCase.pl $pairwise | sort -k2,2n > $prefix.tsv",{jobname=>"eigen",numcpus=>$$settings{numcpus}});
}


######
## phylogeny subroutines
#####


sub phylogenies{
  my($inAln,$settings)=@_;

  my $sge=Schedule::SGELK->new(keep=>1,numcpus=>$$settings{numcpus});
  my $informativeAln=$inAln; # if an informative MSA is not specified, then this one will do.
  $informativeAln=removeUninformativeSites($inAln,$$settings{alnPrefix},$sge,$settings) if($$settings{alnPrefix});
  my $tree=inferPhylogeny($informativeAln,$$settings{treePrefix},$sge,$settings) if($$settings{treePrefix});

  $sge->wrapItUp();  # final wrap-up before leaving this sub
}

sub removeUninformativeSites{
  my($inAln,$outPrefix,$sge,$settings)=@_;
  my $informative="$outPrefix.aln.fas";
  if(-f $informative){
    logmsg "$informative was found.  I will not recalculate.";
    return $informative;
  }
  logmsg "Removing uninformative sites from the alignment and putting it into $informative";
  logmsg "  removeUninformativeSites.pl < '$inAln' > '$informative'";
  $sge->pleaseExecute("removeUninformativeSites.pl < '$inAln' > '$informative'",{jobname=>"removeUninformativeSites",numcpus=>1});
  return $informative;
}

sub inferPhylogeny{
  my($inAln,$prefix,$sge,$settings)=@_;
  my $treeFile="$prefix.RAxML_bipartitions";
  if(-f "$prefix.RAxML_bipartitions"){
    logmsg "$prefix.RAxML_bipartitions was found. I will not recalculate the phylogeny.";
    return $treeFile;
  }
  system("rm -fv $$settings{tempdir}/RAxML*");
  logmsg "Running raxml";
  logmsg "  cd $$settings{tempdir}; launch_raxml.sh $inAln suffix";
  #$inAln=File::Spec->rel2abs($inAln);
  my $rand=int(rand(99999));
  $sge->pleaseExecute("cd $$settings{tempdir}; launch_raxml.sh ../$inAln suffix",{jobname=>"raxml$rand",numcpus=>$$settings{numcpus}});

  # Move those files over when finished
  $sge->wrapItUp();
  for (qw(RAxML_bestTree RAxML_bipartitionsBranchLabels RAxML_bipartitions RAxML_bootstrap RAxML_info)){
    logmsg "I will move $$settings{tempdir}/$_.suffix to $prefix.$_";
    $sge->pleaseExecute("mv -v $$settings{tempdir}/$_.suffix $prefix.$_",{jobname=>"mv$_"});
    #$sge->pleaseExecute("mv -v $$settings{tempdir}/$_.suffix $prefix.$_",{-hold_jid=>"raxml$rand",jobname=>"mv$_"});
  }
  return $treeFile
}

#### 
# things that depend on pairwise and trees
####

sub Fst{
  my($inAln,$pairwisePrefix,$treePrefix,$fstPrefix,$settings)=@_;
  my $fstTree="$fstPrefix.fst.dnd";
  if(-f $fstTree){
    logmsg "$fstTree fst tree was found.  Not recalculating.";
    return $fstTree;
  }
  my $sge=Schedule::SGELK->new(keep=>1,numcpus=>$$settings{numcpus});
  my @command=(
        "applyFstToTree.pl --numcpus $$settings{numcpus} -t $treePrefix.RAxML_bipartitions -p $pairwisePrefix.tsv --outprefix $fstPrefix --outputType averages > $fstPrefix.avg.tsv",
        "applyFstToTree.pl --numcpus $$settings{numcpus} -t $treePrefix.RAxML_bipartitions -p $pairwisePrefix.tsv --outprefix $fstPrefix --outputType samples > $fstPrefix.samples.tsv",
  );
  for (@command){
    logmsg "COMMAND  $_";
    $sge->pleaseExecute($_,{numcpus=>$$settings{numcpus}});
  }
  $sge->wrapItUp();
  return $fstTree;
}

sub usage{
  local $0=basename $0;
  "Process an MSA from a hqSNP pipeline and get useful information
  Usage: $0 file.fasta
  -g   groups.txt       Text files of genome names, one per line (multiple -g are encouraged).
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
