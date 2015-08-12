#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/sum/;
use File::Basename qw/fileparse basename/;
use Bio::Perl;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg @bamExt @fastaExt/;
my $fastaExtRegex=join('$|',@fastaExt).'$';
my $bamExtRegex=join('$|',@bamExt).'$';

local $0=fileparse $0;
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help chunksize=i numchunks=i));
  die usage() if($$settings{help});
  $$settings{chunksize}||=0;
  if(!$$settings{chunksize}){
    $$settings{numchunks}||=32;
  }

  my $infile=$ARGV[0];
  die "ERROR: need infile file!\n".usage() if(!$infile);

  my $contigLength={};
  my($name,$path,$suffix)=fileparse($infile,@bamExt,@fastaExt);
  if($suffix=~/$bamExtRegex/){
    $contigLength=bamLengths($infile,$settings);
  } elsif($suffix=~/$fastaExtRegex/){
    $contigLength=fastaLengths($infile,$settings);
  } else {
    die "ERROR: I do not understand extension $suffix";
  }

  printChunks($contigLength,$settings);

  return 0;
}

sub fastaLengths{
  my($fasta,$settings)=@_;
  my $in=Bio::SeqIO->new(-file=>$fasta,-format=>"fasta");

  my %length;
  while(my $seq=$in->next_seq){
    $length{$seq->id}=$seq->length;
  }
  return \%length;
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
  
  # If the number of chunks is set, then set the chunk length
  # equal to the size divided by the number desired.
  if($$settings{numchunks}){
    my $totalLength=sum(values(%$contigLength));
    $$settings{chunksize}=int($totalLength/$$settings{numchunks})+1;
  }

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
  "Creates a list of regions of a bam file, suitable for samtools' and bcftools' regions flag
  Usage: $0 infile > regions.txt
  infile must either have an extension of either .bam or .fasta

  --chunksize  0   The size of each region.
  --numchunks  32  If chunksize is not set, how many chunks should there be?
  "
}

