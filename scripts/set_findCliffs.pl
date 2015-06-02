#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use List::Util qw/sum max min/;
use Statistics::Descriptive;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use threads;
use Thread::Queue;
use threads::shared;
use POSIX qw/ceil/;
use Statistics::LineFit;

use lib "$FindBin::RealBin/../lib";
use Number::Range;
use LyveSET qw/logmsg/;

my $readBamStick :shared; # A thread must hold the stick before it can read the bam

$0=basename($0);
exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i tempdir=s)) or die $!;
  
  my $bam=$ARGV[0] || "";
  die "ERROR: bam file not given!\n".usage() if(!$bam);
  die "ERROR: could not find $bam\n".usage() if(!-e $bam);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is $$settings{tempdir}";

  my @region=findCliffs($bam,$settings);

  return 0;
}

sub findCliffs{
  my($bam,$settings)=@_;

  # Find the depth of the whole bam file one time, to save cpu.
  # This can be done while regions are being defined in the next loop.
  logmsg "Calculating depth of $bam";
  my $depthfile="$$settings{tempdir}/".basename($bam).".depth";
  system("samtools depth $bam > $depthfile");
  die "ERROR with samtools depth: $!" if $?;
  my %depth;
  open(DEPTH,$depthfile) or die $!;
  while(<DEPTH>){
    chomp;
    my($seqname,$pos,$d)=split /\t/;
    $depth{$seqname}{$pos}=$d;
  }
  close DEPTH;

  # Generate the regions that will be looked at.
  my $windowSize=250;
  my $stepSize=int($windowSize/2);
  my $seqinfo=getSeqInfo($bam,$settings);
  my @region; # set of coordinates to analyze for depth
  my %region; # hash of different range objects for cliffs
  while(my($seqname,$length)=each(%$seqinfo)){
    $region{$seqname}=Number::Range->new;
    for(my $pos=1;$pos<=$length;$pos+=$stepSize){
      my $stop=$pos+$windowSize-1;
      push(@region,[$bam,$seqname,$pos,$stop]);
    }
  }
  # Split the windows into thread-sized chunks
  logmsg "Defining regions to give to the cliff-detector";
  my $numPerThread=ceil(@region/$$settings{numcpus});
  my @regionSet;
  my @depthSet;
  for(my $i=0;$i<$$settings{numcpus};$i++){
    # Split up the regions
    $regionSet[$i]=[splice(@region,0,$numPerThread)];

    # Also split up the depth hashes
    for(@{$regionSet[$i]}){ # [bam,seqname,start,stop]
      my(undef,$seqname,$start,$stop)=@$_;
      $depthSet[$i]{$seqname}{$_}=$depth{$seqname}{$_}||0 for($start..$stop);
    }
  }
  logmsg "Now investigating those regions for cliffs.";

  # Make threads and distribute the regions array evenly to them all
  my @thr;
  for my $i(0..$$settings{numcpus}-1){
    my %depthCopy=%{$depthSet[$i]};
    $thr[$i]=threads->new(\&findCliffsInRegionWorker,$regionSet[$i],\%depthCopy,$settings);
  }

  # Gather results from the threads by making
  # a set of ranges
  my @cliff;
  for my $t(@thr){
    logmsg "Combining ranges found in TID".$t->tid;
    my $cliffs=$t->join;
    push(@cliff,@$cliffs);
  }

  # Combine ranges
  my %cliffRange;
  for my $cliff(@cliff){
    my $name=$$cliff{region}[0];
    $cliffRange{$name}||=Number::Range->new;
    no warnings;
    $cliffRange{$name}->addrange($$cliff{region}[1]..$$cliff{region}[2]);
  }

  # Print ranges
  while(my($seqname,$cliff)=each(%cliffRange)){
    my $range=$cliff->range();

    # some internal sorting by start
    my @pos;
    while($range=~/(\d+)\.\.(\d+),?/g){
      push(@pos,[$1,$2]);
    }

    for my $lohi(sort {$$a[0] <=> $$b[0]} @pos){
      my($lo,$hi)=@$lohi;
      my $name="$seqname:$lo";
      print join("\t",$seqname,$lo,$hi,$name)."\n";
    }
  }

}

sub findCliffsInRegionWorker{
  my($regionArr,$depth,$settings)=@_;

  logmsg "Init";

  my @region; # array of hashes
  my $i=0;
  for my $tmp(@$regionArr){
    my($bam,$seqname,$pos,$lastPos)=@$tmp;
    my $r=findCliffsBySlope($bam,$depth,$seqname,$pos,$lastPos,$settings);
    push(@region,@$r);
    $i++;
  }
  logmsg "Finished looking at $i regions";
  return \@region;
}

sub readBamDepth{
  my($bam,$region,$settings)=@_;
  lock $readBamStick;
  open(DEPTH,"samtools depth -r '$region' $bam | ") or die "ERROR: could not open $bam:$!";
  my %depth;
  while(<DEPTH>){
    chomp;
    my($seqname,$pos,$depth)=split /\t/;
    $depth{$seqname}{$pos}=$depth;
  }
  close DEPTH;

  return \%depth;
}

sub findCliffsBySlope{
  my($bam,$allDepth,$seqname,$pos,$lastPos,$settings)=@_;
  
  $$settings{slopeThreshold}||=3; # detect any region that has a slope bigger than this abs number
  
  my $windowSize=15;
  my $stepSize=int($windowSize/2);

  # Figure out the slope every 15 bp in this region from $pos to $lastPos.
  # If it goes beyond the threshold then it is a cliff.
  my @cliff=();
  while($pos<=$lastPos){

    # Load up information for line-fitting (ie, trendline)
    my (@depth,@pos); # y and x on the linefit graph
    my $start=$pos;
    my $stop=$start+$windowSize-1;
    while($pos<=$stop){
      
      my $depth=$$allDepth{$seqname}{$pos} || 0;
      push(@depth,$depth);
      push(@pos,$pos);

      $pos++;
    }

    # Get the actual linefit statistics
    my $linefit=Statistics::LineFit->new();
    $linefit->setData(\@pos,\@depth) or die "Invalid depths to figure out the slope! \n @depth\n";

    # Don't worry about this being a cliff if rSquared is not significant
    next if(!defined($linefit->rSquared()) || $linefit->rSquared() < 0.5);

    # Get the slope of the line
    my($intercept,$slope)=$linefit->coefficients();

    # Record this as a cliff if the slope is steep enough
    next if(abs($slope) < $$settings{slopeThreshold});
    push(@cliff,{
      region=>[$seqname,$start,$stop],
      intercept=>$intercept,
      slope=>$slope,
      hi=>max(@depth),
      lo=>min(@depth),
      name=>"$seqname:$start",
      rSquared=>$linefit->rSquared(),
    });


  }

  return \@cliff;
}

# Finds sequence names (contigs) and their lengths
sub getSeqInfo{
  my($bam,$settings)=@_;
  my (%seqinfo);
  open(BAMHEADER,"samtools view -H $bam | ") or die "ERROR: could not open $bam with samtools view -H: $!";
  while(<BAMHEADER>){
    chomp;
    next if(!/\@SQ/);
    my($SQ,@fields)=split(/\t/);
    my %F;
    for my $keyvalue(@fields){
      my($key,$value)=split(/:/,$keyvalue);
      $F{$key}=$value;
    }
    $seqinfo{$F{SN}}=$F{LN} if($F{SN} && $F{LN});
  }
  close BAMHEADER;
  return \%seqinfo;
}

sub depth{
  my($bam,$region,$settings)=@_;
  open(DEPTH,"samtools depth -r '$region' $bam | ") or die "ERROR: could not open $bam with samtools depth: $!";
  my @depth;
  while(<DEPTH>){
    chomp;
    my($seqname,$pos,$depth)=split /\t/;
    push(@depth,$depth);
  }
  close DEPTH;
  
  return 0 if(!@depth);
  return sum(@depth)/scalar(@depth);
}

sub usage{
  "Finds cliffs in the depth in a given sorted/indexed bam file.
  Usage: $0 file.bam
  --numcpus 1
  "
}
