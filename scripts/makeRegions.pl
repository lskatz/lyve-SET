#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/sum/;
use File::Basename qw/fileparse basename/;
use Bio::Perl;

use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg @bamExt @fastaExt @vcfExt/;
my $fastaExtRegex=join('$|',@fastaExt).'$';
my $bamExtRegex=join('$|',@bamExt).'$';
my $vcfExtRegex=join('$|',@vcfExt).'$';

local $0=fileparse $0;
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help chunksize=i numchunks=i numcpus=i minlength=i overlapby=i));
  die usage() if($$settings{help});
  $$settings{minlength}||=1;
  $$settings{chunksize}||=0;
  $$settings{overlapby}||=0;
  if(!$$settings{chunksize}){
    $$settings{numchunks}||=1;
  }
  $$settings{numcpus}||=1;

  die "ERROR: need infile file!\n".usage() if(!@ARGV);

  # Start off multithreading
  my $fileQ=Thread::Queue->new(@ARGV);
  my @thr;
  for my $i(0..$$settings{numcpus}-1){
    $thr[$i]=threads->new(\&contigLengthsWorker,$fileQ,$settings);
  }

  for my $infile(@ARGV){
    my $contigLength={};
    my($name,$path,$suffix)=fileparse($infile,@bamExt,@fastaExt,@vcfExt);
    logmsg "Enqueuing $infile with extension $suffix";
    $fileQ->enqueue($infile);
  }
  $fileQ->enqueue(undef) for(@thr);

  # Combine contig length results from threads
  my %cLength;
  for(@thr){
    my $tmp=$_->join;
    %cLength=(%cLength,%$tmp);
  }

  # print combined results
  printChunks(\%cLength,$settings);

  return 0;
}

sub contigLengthsWorker{
  my($fileQ,$settings)=@_;
  my $cLength;
  while(defined(my $infile=$fileQ->dequeue)){
    my $contigLength={};
    my($name,$path,$suffix)=fileparse($infile,@bamExt,@fastaExt,@vcfExt);
    if($suffix=~/$bamExtRegex/){
      $contigLength=bamLengths($infile,$settings);
    } elsif($suffix=~/$fastaExtRegex/){
      $contigLength=fastaLengths($infile,$settings);
    } elsif($suffix=~/$vcfExtRegex/){
      $contigLength=vcfLengths($infile,$settings);
    } else {
      die "ERROR: I do not understand extension $suffix from filename $infile";
    }
    
    # Combine
    while(my($seqname,$length)=each(%$contigLength)){
      $$cLength{$seqname}=$length if(!defined($$cLength{$seqname}));
      $$cLength{$seqname}=$length if($$cLength{$seqname} < $length);
    }

  }
  return $cLength;
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

# Finding the length of a VCF is sort of difficult.
# There are three different ways to find it:
#   1. The ##reference line, and then parsing the fasta file
#   2. The bcftools index -s function
#   3. Parsing the file to find the max coordinate per seqname (slow)
sub vcfLengths{
  my($vcf,$settings)=@_;
  
  # See if ##reference=reference can be found, and if so,
  # if this script can parse that fasta file for lengths.
  my $fasta;
  open(VCFGZ,"zgrep '^##reference' $vcf | ") or die "ERROR: could not open $vcf for reading: $!";
  while(<VCFGZ>){
    if(/^##reference=(.+)/){
      $fasta=$1;
      last if(-e $fasta);

      logmsg "Found reference genome $fasta in the vcf file but it does not exist where I expect it";
      $fasta="";
    }
  }
  close VCFGZ;
  if($fasta){
    return fastaLengths($fasta,$settings);
  }

  my %max; # used for any of the next methods

  # See if bcftools can help find the length if the ref is not found
  my $bcftoolsIndex=`bcftools index -s '$vcf' 2>/dev/null`;
  #logmsg "ERROR with command => bcftools index -s '$vcf'\n  $0 @ARGV" if $?;
  for(split(/\n/,$bcftoolsIndex)){
    chomp;
    my($seqname,$start,$stop)=split /\t/;
    $max{$seqname}=$stop;
  }
  close VCFGZ;
  return \%max if(keys(%max)>0);

  #die "ERROR: could not determine contig lengths for $vcf using either the ##reference tag nor 'bcftools index -s'";

  #die Dumper "ERROR",\%max;

  # If a reference is not found and bcftools didn't come through, 
  # then just find the coordinates in the file itself (slow step).
  # Uses a sort step such that the highest coordinate of every 
  # contig is listed first for each contig.
  open(VCFGZ,"zcat $vcf | cut -f 1,2 | sort -k1,1 -k2,2nr |") or die "ERROR: could not open $vcf for reading: $!";
  while(<VCFGZ>){
    next if(/^#/);
    my($seqname,$pos)=split /\t/;
    next if(defined($max{$seqname}));
    $max{$seqname}=$pos+0;
  }
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

  # Print out the regions
  my $chunksize=$$settings{chunksize};
  while(my($seqname,$length)=each(%$contigLength)){
    for(my $start=1;$start<$length;$start+=$chunksize){
      my $end = $start + $chunksize - 1 + $$settings{overlapby};
      $end=$length if($end > $length);
      print "$seqname:$start-$end\n";
    }
  }

}

sub usage{
  "Creates a list of regions of a bam file, suitable for samtools' and bcftools' regions flag
  Usage: $0 infile [infile2...] > regions.txt
  infile must either have an extension of either .bam, .vcf.gz, or .fasta.  Bam and vcf files must be indexed.
  It's a good idea to give multiple vcf files since they don't really betray where the last coordinate is. In other words, another vcf file might have a snp beyond where the curren vcf's snps are

  --chunksize  0  The size of each region.
  --numchunks  1  If chunksize is not set, how many chunks should there be?
                  NOTE: Despite what is requested, there will be at least 
                  one chunk per contig.
  --overlapby  0  How many bp each region should overlap by
  --minlength  1  The minimum number of base pairs allowed per chunk
  "
}

