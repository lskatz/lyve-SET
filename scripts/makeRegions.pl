#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Bio::Perl;

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help chunksize=i));
  die usage() if($$settings{help});
  $$settings{chunksize}||=10000;

  my $bam=$ARGV[0];
  die "ERROR: need bam file!\n".usage() if(!$bam);

  my $contigLength=lengths($bam,$settings);
  printChunks($contigLength,$settings);

  return 0;
}

sub lengths{
  my($sorted,$settings)=@_;
  # Find the total length of the reference genome
  logmsg "Finding the total length of each contig";
  my %max;
  open(BAMHEADER,"samtools view -H $sorted | ") or die "ERROR: could not use samtools view on $sorted: $!";
  while(<BAMHEADER>){
    chomp;
    my($type,@F)=split /\t/;
    $type=~s/^\@//; # remove header indicator
    next if($type ne 'SQ');
    my %F;
    for(@F){
      my($key,$value)=split(/:/);
      $F{$key}=$value;
    }
    $max{$F{SN}}=$F{LN};
  }
  close BAMHEADER;
  return \%max;
}

sub printChunks{
  my($contigLength,$settings)=@_;
  
  my $chunksize=$$settings{chunksize};
  while(my($seqname,$length)=each(%$contigLength)){
    for(my $start=1;$start<$length;$start+=$chunksize){
      my $end = $start + $chunksize - 1;
      $end=$length if($end > $length);
      print "$seqname:$start-$end\n";
    }
  }

}

sub usage{
  "Creates a list of regions of a bam file, suitable for samtools and bcftools
  Usage: $0 file.sorted.bam > regions.txt
  --chunksize  10000  The size of each region
  "
}

