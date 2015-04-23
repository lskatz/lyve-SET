#!/usr/bin/env perl
# Finds the kmer jaccard distance between any two sets of raw reads
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Math::BigInt;

sub logmsg{ $|++; print STDERR "@_\n"; $|--;}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help quiet coverage=i kmerlength=i numcpus=i)) or die $!;
  die usage() if(!@ARGV || $$settings{help});
  my @asm=@ARGV;
  $$settings{coverage}||=2; 
  $$settings{coverage}=1 if($$settings{coverage}<1);
  $$settings{kmerlength}||=18;


  jaccardDistance(\@asm,$settings);

  return 0;
}

sub jaccardDistance{
  my($genome,$settings)=@_;
  my %jDist;
  for(my $i=0;$i<@$genome-1;$i++){
    my %k1=kmerCount($$genome[$i],$settings);
    for(my $j=$i+1;$j<@$genome;$j++){
      my %k2=kmerCount($$genome[$j],$settings);
      my $jDist=jDist(\%k1,\%k2,$settings);
      print join("\t",@$genome[$i,$j],$jDist)."\n";
    }
  }
}

sub jDist{
  my($k1,$k2,$settings)=@_;
  my $minKCoverage=$$settings{coverage};

  logmsg "Finding intersection and union of kmers";
  my %kmerSet=kmerSets($k1,$k2,$settings);

  my $jDist=1-($kmerSet{intersection} / $kmerSet{union});

  logmsg "$jDist=1-($kmerSet{intersection} / $kmerSet{union})";
  return $jDist;
}

sub kmerSets{
  my($k1,$k2,$settings)=@_;
  
  my($intersectionCount,%union);

  # Find uniq kmers in the first set of kmers.
  # Also find the union.
  for my $kmer(keys(%$k1)){
    $intersectionCount++ if(!$$k2{$kmer});
    $union{$kmer}=1;
  }

  # Find uniq kmers in the second set of kmers.
  # Also find the union.
  for my $kmer(keys(%$k2)){
    $intersectionCount++ if(!$$k1{$kmer});
    $union{$kmer}=1;
  }

  my $unionCount=scalar(keys(%union));

  return (intersection=>$intersectionCount,union=>$unionCount);
}

sub kmerCount{
  my($genome,$settings)=@_;
  my $kmerLength=$$settings{kmerlength};
  my $minKCoverage=$$settings{coverage};
  logmsg "Counting $kmerLength-mers for $genome";
  my($name,$path,$suffix)=fileparse($genome,qw(.fastq.gz .fastq));
  if($suffix=~/\.fastq\.gz$/){
    open(FILE,"gunzip -c '$genome' |") or die "I could not open $genome with gzip: $!";
  } elsif($suffix=~/\.fastq$/){
    open(FILE,"<",$genome) or die "I could not open $genome: $!";
  } else {
    die "I do not understand the extension on $genome";
  }

  # count kmers
  my %kmer;
  my $i=0;
  while(my $read=<FILE>){
    $i++;
    logmsg "Finished with ".($i/4)." reads" if($i % 1000000 == 0);
    next if($i % 4 != 2);
    chomp $read;

    # How long to read along the sequence
    my $length=length($read)-$kmerLength+1;
    # Start saving on memory by converting the read immediately to a number.
    $read=~tr/ATCGNatcgn/0123401234/;
    my @readKmer=(); 
    for(my $j=0;$j<$length;$j++){
      # Mess around with this smaller kmer array per read
      # and hopefully save on cpu/memory.
      push(@readKmer,
        "1".substr($read,$j,$kmerLength) # make sure there is no zero in front
      );
    }

    # Now put back all those kmers into the large hash.
    for my $kmer(@readKmer){
      #$kmer{Math::BigInt->new($kmer)}++;
      $kmer=$kmer+0; # convert to int
      $kmer{$kmer}++;
    }
  }
  close FILE;

  # remove kmers with low depth
  while(my($kmer,$count)=each(%kmer)){
    delete($kmer{$kmer}) if($count<$minKCoverage);
  }

  logmsg "Found ".scalar(keys(%kmer))." unique kmers of depth > $minKCoverage";
  return %kmer;
}

sub usage{
  "Finds the jaccard distance between two assemblies using mummer. With more genomes, it creates a table.
  Jaccard distance: (kmer method) counts 18-mers and calculates 1 - (intersection/union)
  Usage: $0 file.fastq.gz file2.fastq.gz [file3.fastq.gz ...]
  -q for minimal stdout
  -c minimum kmer coverage before it counts. Default: 2
  -k kmer length. Default: 18
  -n numcpus
  "
}
