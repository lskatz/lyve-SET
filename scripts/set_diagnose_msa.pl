#!/usr/bin/env perl
# Diagnoses potential problems in a SET multiple sequence alignment
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Perl;
use Bio::AlignIO;
use File::Basename qw/basename/;
use POSIX qw/ceil/;

use FindBin;
use lib "$FindBin::RealBin/../lib";
#use Statistics::Normality qw/shapiro_wilk_test dagostino_k_square_test/;
use Statistics::Descriptive;

local $0=basename $0;
sub logmsg {
  my $FH = *STDOUT;
  my $sub=(caller(1))[3];
  $sub=~s/main:://;
  print $FH "$0:$sub: @_\n";
}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help ambiguity-tolerance=i msa=s)) or die $!;
  $$settings{'ambiguity-tolerance'}||=10;
  $$settings{msa}||="out.aln.fas";
  die usage() if($$settings{help});

  my($dir)=@ARGV;
  $dir||=".";

  if(!-e $dir){
    logmsg "ERROR: Could not find the MSA directory $dir\n".usage();
  }

  my $msa="$dir/$$settings{msa}";
  die "ERROR: Could not find MSA $msa\n".usage() if(!-e $msa);
  my $msaObj=Bio::AlignIO->new(-file=>$msa)->next_aln;

  findAmbiguities($msaObj,$settings);
  numberOfSitesMasked($msaObj,$settings);

  return 0;
}

sub findAmbiguities{
  my($msaObj,$settings)=@_;
  my $length=$msaObj->length;
  for my $seq($msaObj->each_seq){
    my $sequence=$seq->seq;
    my $N=($sequence=~s/[^ATCG]//gi);
    $N||=0;
    my $percentAmbiguous=$N/$length * 100;

    if($percentAmbiguous > $$settings{'ambiguity-tolerance'}){
      $percentAmbiguous=sprintf("%0.1f",$percentAmbiguous);
      logmsg "WARNING: ".$seq->id." has $percentAmbiguous% ambiguous nucleotides";
    }
  }
  return;
}

sub numberOfSitesMasked{
  my($msaObj,$settings)=@_;
  # TODO make this a parameter with a better name
  $$settings{frequencyMaskedDefinesMaskedSite}||=0.5;
  my $numSeqs=$msaObj->num_sequences;
  my $numSeqsMeansMasked=ceil($numSeqs * $$settings{frequencyMaskedDefinesMaskedSite});

  # Count how many genomes per site are masked
  my $length=$msaObj->length;
  my @maskedSitesCounter;
  for my $seq($msaObj->each_seq){
    # Remember the DNA sequence to speed things up
    my $sequence=$seq->seq;
    for(my $i=0;$i<$length;$i++){
      if(substr($sequence,$i,1) =~/[^atcg]/i){
        $maskedSitesCounter[$i]++;
      }
    }
  }

  # Find the number of sites that are masked above a 
  # threshold percentage of genomes.
  my $numMaskedSites=0;
  my $numHasAnyMasking=0;
  my $numTotallyMasked=0;
  for (@maskedSitesCounter){
    $_||=0;
    if($_ >= $numSeqsMeansMasked){
      $numMaskedSites++;
    } 
    if($_ > 0){
      $numHasAnyMasking++;
    }
    if($_ == $numSeqs){
      $numTotallyMasked++;
    }
  }

  # Make a report
  my $percentMasked=int($numMaskedSites/$length * 100);
  logmsg "The MSA is $percentMasked% masked, defined by sites that are masked with a frequency > $$settings{frequencyMaskedDefinesMaskedSite}";
  my $percentAtAllMasked=int($numHasAnyMasking/$length * 100);
  logmsg "$percentAtAllMasked% of the sites in the MSA have at least one masked nucleotide";
  my $percentTotallyMasked=int($numTotallyMasked/$length * 100);
  logmsg "$percentTotallyMasked% of the sites in the MSA are totally masked and are unnecessary for the analysis";

  return;
}

sub outliers{
  #print( (ref($_)||"INT")."\n") for(@_);
  my($stat,$datum)=@_;
  ...;
}

sub usage{
  "Diagnoses any potential problems with the SET multiple sequence alignment directory
  Usage: $0 [project/msa]
  project/msa  The MSA directory. Default: current directory
  --msa                  out.aln.fas
  --ambiguity-tolerance  10
  --help                 Help menu
  "
}
