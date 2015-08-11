#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/fileparse basename/;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg @fastaExt/;

exit(main());

sub main{
  my $settings={};
  die usage($settings) if(!@ARGV);
  GetOptions($settings,qw(help reference=s project=s matrix=s));
  $$settings{maskedThresholdPercent}=10;
  $$settings{project}||="";
  $$settings{reference}||="";
  $$settings{matrix}||="";

  die usage($settings) if($$settings{help});

  # check that this is a Lyve-SET project if one is given
  if($$settings{project}){
    # See if this project exists and is legit
    system("set_manage.pl $$settings{project}");
    die "ERROR: $$settings{project} is not a Lyve-SET project!" if $?;

    # Now find the reference genome
    # Make an array of possible reference genomes
    my @fasta=("$$settings{project}/reference/reference.fasta");
    for my $ext(@fastaExt){
      push(@fasta,(glob("$$settings{project}/reference/*.$ext")));
    }
    # Find the first one that exists and has a nonzero filesize
    for my $fasta(@fasta){
      if(-e $fasta && -s $fasta > 0){
        $$settings{reference} ||= $fasta;
        logmsg "Found $$settings{reference} reference genome";
        last;
      }
    }

    # Find the SNP matrix
    my $msadir="$$settings{project}/msa";
    my @tsv=("$msadir/out.snpmatrix.tsv","$msadir/out.bcftoolsquery.tsv","$msadir/out.filteredMatrix.tsv");
    for my $tsv(@tsv){
      if(-e $tsv && -s $tsv > 0){
        $$settings{matrix}||=$tsv;
        logmsg "Found $$settings{matrix} SNP matrix";
        last;
      }
    }
    
  }
  # END if project is given

  # If no project is given, then assign based on what the user specified
  else {
    ($$settings{matrix})=@ARGV;
  }

  # Now do some diagnostics
  my $refInfo=readAssembly($settings);
  reportMaskedGenomes($$settings{matrix},$refInfo,$settings);

  return 0;
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
      die "ERROR: field $_ was not found in $matrix" if(!$site{$_});
      $site{$_}=$sample{$_};
      delete($sample{$_});
    }
    my $REF=$site{REF}; # needed for substitution in the loop below
    
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
  my $numHqSites      =$numSites - $masked100;
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
    my $percentHq=int($numHqSites/$asmLength*10000)/100;
    logmsg "With $minLengthStr min length contigs, $percentHq% is represented in the matrix";
  }
}

sub usage{
  my($settings)=@_;
  local $0=basename $0;
  my $help = "Diagnoses a SNP matrix file
  USAGE: $0 -p set-project
  note: set-project must be a Lyve-SET project folder
  -h for more help
  ";
  return $help if(!$$settings{help});
  $help.="
  Alternate USAGE: $0 --matrix snpmatrix.tsv [-r reference.fasta]
  note: snpmatrix.tsv must have as the first three columns: CHROM, POS, REF and must have subsequent columns pertaining to each sample
  ";
  return $help;
}
