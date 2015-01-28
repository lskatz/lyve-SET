#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Bio::Perl;

use FindBin;
# Include samtools and bcftools
$ENV{PATH}=$ENV{PATH}.":$FindBin::RealBin/../lib/samtools-1.1:$FindBin::RealBin/../lib/samtools-1.1/misc";
$ENV{PATH}=$ENV{PATH}.":$FindBin::RealBin/../lib/bcftools-1.1";


local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={clean=>0};
  GetOptions($settings,qw(help));
  die usage() if($$settings{help});

  my $bam=$ARGV[0];
  die "ERROR: need bam file!\n".usage() if(!$bam);

  bamDepth($bam,$settings);

  return 0;
}

sub bamDepth{
  my($sorted,$settings)=@_;

  # samtools depth output; store it in a hash
  logmsg "Getting coverage levels at each position";
  my %depth;
  open(BAMDEPTH,"samtools depth $sorted | ") or die "ERROR: could not open the bam file with samtools: $!";
  while(<BAMDEPTH>){
    chomp;
    my @F=split /\t/;
    $depth{$F[0]}{$F[1]}=$F[2];
  }
  close BAMDEPTH;
 
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

  # get zero-depth information in there too. First find the total length of each contig.
  logmsg "Backfilling for zero-depth positions";
  open(FIXEDDEPTH,">","$sorted.depth") or die "ERROR: could not write to $sorted.depth: $!";
  for my $contig(keys(%max)){
    for my $i(1..$max{$contig}){
      $depth{$contig}{$i}||=0; 
      print FIXEDDEPTH join("\t",$contig,$i,$depth{$contig}{$i})."\n";
    }
  }
  close FIXEDDEPTH;

  # compress the depth
  logmsg "Compressing";
  system("gzip -v9 $sorted.depth");
  die "ERROR with gzip" if $?;
  
  return 1;
}

sub usage{
  "Runs samtools depth and fills in missing zeros
  Usage: $0 file.sorted.bam 
    Output file: file.sorted.bam.depth[.gz]
  "
}

