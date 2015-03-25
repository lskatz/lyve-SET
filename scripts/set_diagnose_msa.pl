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

use FindBin;
use lib "$FindBin::RealBin/../lib";
use Statistics::Normality qw/shapiro_wilk_test dagostino_k_square_test/;
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
