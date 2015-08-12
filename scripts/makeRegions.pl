#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/fileparse basename/;
use Bio::Perl;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg @bamExt @fastaExt/;
my $fastaExtRegex=join('$|',@fastaExt).'$';
my $bamExtRegex=join('$|',@bamExt).'$';

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help chunksize=i));
  die usage() if($$settings{help});
  $$settings{chunksize}||=10000;

  my $infile=$ARGV[0];
  die "ERROR: need infile file!\n".usage() if(!$infile);

  my $contigLength={};
  my($name,$path,$suffix)=fileparse($infile,@bamExt,@fastaExt);
  if($suffix=~/$bamExtRegex/){
    $contigLength=bamLengths($infile,$settings);
  } elsif($suffix=~/$fastaExtRegex/){
    ...;
  } else {
    die "ERROR: I do not understand extension $suffix";
  }
  printChunks($contigLength,$settings);

  return 0;
}

sub bamLengths{
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
  Usage: $0 infile > regions.txt
  infile must either have an extension of either .bam or .fasta

  --chunksize  10000  The size of each region
  "
}

