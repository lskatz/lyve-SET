#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use List::Util qw/sum/;
use Statistics::Descriptive;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use threads;
use Thread::Queue;
use threads::shared;
use POSIX qw/ceil/;

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
  logmsg "Calculating depth of $bam in the background.";
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

  # Generate the regions that will be looked at
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
  logmsg "Done defining regions.";

  # Make threads and distribute the regions array evenly to them all
  my @thr;
  for my $i(0..$$settings{numcpus}-1){
    my %depthCopy=%{$depthSet[$i]};
    $thr[$i]=threads->new(\&findCliffsInRegionWorker,$regionSet[$i],\%depthCopy,$settings);
  }

  # Gather results from the threads by making
  # a set of ranges
  for my $t(@thr){
    logmsg "Combining ranges found in TID".$t->tid;
    my $ranges=$t->join;
    # Add whatever coordinates the analysis showed us
    for my $setOfRanges(@$ranges){
      next if(!@$setOfRanges);
      for my $c(@$setOfRanges){
        no warnings; # b/c $_->addrange will produce a warning if there are repeated coordinates
        $region{$$c{region}[0]}->addrange($$c{region}[1]..$$c{region}[2]);
      }
    }
  }
  # DONE gathering and combining ranges

  # Print the results in BED format
  for my $seqname(keys(%$seqinfo)){
    next if(!$region{$seqname});
    my $rangeStr=$region{$seqname}->range;
    while($rangeStr=~/(\d+)\.\.(\d+),?/g){
      my($start,$stop)=($1,$2);
      my $name=join(":","cliff",$seqname,"$start-$stop");
      print join("\t",$seqname,$start,$stop,$name)."\n";
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
    my @r=findCliffsInRegion($bam,$depth,$seqname,$pos,$lastPos,$settings);
    push(@region,@r);
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

sub findCliffsInRegion{
  my($bam,$allDepth,$seqname,$pos,$lastPos,$settings)=@_;
  my %region; # $region{contig}=>[ 
              #                    [start,stop],
              #                    [start,stop],
              #                  ]

  my $windowSize=15;
  my $stepSize=int($windowSize/2);

  # Find the average of average depths
  #my $allDepth=readBamDepth($bam,"$seqname:$pos-$lastPos",$settings);
  my %depth;
  for(my $start=$pos;$start<$lastPos;$start+=$stepSize){
    my $stop=$start+$windowSize-1;
    $stop=$lastPos if($stop>$lastPos);
    
    # Calculate this window's depth
    my $total;
    for($start .. $stop){
      $total+=$$allDepth{$seqname}{$_} || 0;
    }
    $depth{$start}=$total/($stop - $start + 1);
    #my $windowDepth=depth($bam,"$seqname:$start-$stop",$settings);
    #$depth{$start}=$windowDepth;
  }
  my $stat=Statistics::Descriptive::Full->new;
  $stat->add_data(values(%depth));
  my($avg,$stdev)=($stat->mean,$stat->standard_deviation);

  # Is anything lower than a cutoff for what I'd expect?
  # Cliffs don't happen all to often, so maybe it's in the bottom 1%
  my $lo=$avg - 3*$stdev;
  for($lo){
    $_=0 if($_<0);
    $_=sprintf("%0.2f",$_);
  }

  # Figure out the cliffs using $lo as a cutoff
  my @region;
  # I don't want to look at the edges of the larger window, and
  # anything I miss will be encompassed in the next cascading
  # window anyway.
  my @sortedStart=sort {$a<=>$b} keys(%depth);
  shift(@sortedStart); pop(@sortedStart);
  for(my $i=0;$i<@sortedStart;$i++){
    my $start=$sortedStart[$i];
    my $depth=$depth{$start};
    next if($depth > $lo); # don't call this a cliff it's not too low
    push(@region,{depth=>$depth,CI=>[$lo,9999],region=>[$seqname,$start,$start+$windowSize]});
  }
  #logmsg "TID:".threads->tid." $seqname:$pos-$lastPos  Depth 95% CI interval: [$lo,$hi]";
  return \@region;
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
