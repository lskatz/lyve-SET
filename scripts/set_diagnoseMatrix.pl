#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help reference=s));
  $$settings{maskedThresholdPercent}=10;

  die usage() if($$settings{help} || !@ARGV);

  my $refInfo=readAssembly($settings);
  for my $matrix(@ARGV){
    reportMaskedGenomes($matrix,$refInfo,$settings);
  }
}

sub readAssembly{
  my($settings)=@_;
  my $asmLength={};
  my $ref=$$settings{reference} || return $asmLength;

  # run assembly metrics on min=150, 250, and 1000 to reflect different meanings
  for my $min(qw(150 250 1000)){
    # run_assembly_metrics.pl -m 150 ../reference/reference.fasta -s genomeLength  --number
    $$asmLength{$min} = `run_assembly_metrics.pl -m $min -s genomeLength --number $ref`;
    die "ERROR: could not find the assembly length in $ref" if $?;
    chomp($$asmLength{$min});
    
  }
  return $asmLength;
}

sub reportMaskedGenomes{
  my($matrix,$refInfo,$settings)=@_;
  open(MATRIX,$matrix) or die "ERROR: could not open $matrix for reading: $!";
  my $header=<MATRIX>;
  $header=~s/^#//g;                # remove leading pound sign
  $header=~s/^\s+|\s+$//g;         # remove whitespace 
  my @header=split(/\t/,$header);
  $_=~s/^\[\d+\]// for(@header);   # remove leading col counter
  $_=~s/:GT// for(@header);        # remove following :GT notation

  # Read the SNP matrix file
  my %maskCounter=();
  my $numSites=0;
  my($maskedSitesCounter, $masked50, $masked100)=(0,0,0);
  while(<MATRIX>){
    chomp;
    my(%sample,%site);

    # index the columns
    my @F=split(/\t/);
    @sample{@header}=@F;

    # split that index into site info and sample info
    for(qw(CHROM POS REF)){
      $site{$_}=$sample{$_};
      delete($sample{$_});
    }
    my $REF=$site{REF};
    
    my $numSamples=scalar(keys(%sample)); # TODO get this outside of the loop
    my $numSamplesMaskedHere=0;
    while(my($sample,$nt)=each(%sample)){
      $nt=~s/\./$REF/;
      # Count masked bases
      if($nt=~/N/i){
        $maskCounter{$sample}++;
        $numSamplesMaskedHere++;
      }
    }

    $maskedSitesCounter++ if($numSamplesMaskedHere>0);
    $masked50++ if($numSamplesMaskedHere/$numSamples >= 0.5);
    $masked100++ if($numSamplesMaskedHere/$numSamples == 1);

    # Keep a running tally of sites
    $numSites++;
  }
  close MATRIX;

  # Percentage of masked bases per genome
  my @genome=sort { $a cmp $b } keys(%maskCounter); # sort, for prettier output
  for my $sample(@genome){
    my $count=$maskCounter{$sample};
    my $percentGenomeMasked=int($count/$numSites*10000)/100;
    logmsg "$sample is $percentGenomeMasked% masked" if($percentGenomeMasked > 10);
  }

  # How many sites are masked
  my $percentMasked   =int($maskedSitesCounter/$numSites*10000)/100;
  my $percentMasked50 =int($masked50/$numSites*10000)/100;
  my $percentMasked100=int($masked100/$numSites*10000)/100;
  logmsg "Sites masked at all: $percentMasked%";
  logmsg "Sites where 50% of the site is masked: $percentMasked50%";
  logmsg "Sites where 100% of the site is masked: $percentMasked100%";

  # What percent of the genome is represented in the matrix?
  my @minLength=sort {$a <=> $b} keys(%$refInfo);
  for my $minLength(@minLength){
  #while(my($minLength,$asmLength)=each(%$refInfo)){
    my $asmLength=$$refInfo{$minLength};
    my $minLengthStr=sprintf("%0.2f%s",int($minLength*100)/100000,"k");
    my $percentHq=int($numSites/$asmLength*10000)/100;
    logmsg "With $minLengthStr min length contigs, $percentHq% is represented in the matrix";
  }
}

sub usage{
  "$0: diagnoses a SNP matrix file
  Usage: $0 snpmatrix.tsv
  NOTE: snpmatrix.tsv must have as the first three columns: CHROM, POS, REF and must have subsequent columns pertaining to each sample
  -r reference.fasta (optional)
  "
}
