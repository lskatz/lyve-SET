#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use List::Util qw/sum max min/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir tempfile/;
use threads;
use Thread::Queue;
use Bio::FeatureIO;

use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;
use lib "$FindBin::RealBin/../lib/lib/perl5";
#use Number::Range;
use Array::IntSpan;
use Statistics::LineFit;
use Statistics::Descriptive;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i tempdir=s slopeThreshold=s rSquaredThreshold=s windowSize=i)) or die $!;

  # Check on parameters
  my $bam=$ARGV[0] || "";
  die "ERROR: bam file not given!\n".usage() if(!$bam);
  die "ERROR: could not find $bam\n".usage() if(!-e $bam);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{slopeThreshold}||=5;      # detect any region that has a slope bigger than this abs number
  $$settings{rSquaredThreshold}||=0.7; # level of significance in the trendline slope
  $$settings{windowSize}||=10;

  logmsg "Temporary directory is $$settings{tempdir}";
  
  # Kick off all the threads.
  # Items being enqueued later go into this subroutine.
  my @thr;
  my $Q=Thread::Queue->new;
  $thr[$_]=threads->new(sub{
    my($Q,$settings)=@_;
    my @cliff;
    while(defined(my $args=$Q->dequeue)){
      my $cliff=findCliffsBySlope(@$args,$settings);
      push(@cliff,@$cliff);
    }
    return \@cliff;
  },$Q,$settings) for(0..$$settings{numcpus}-1);

  # Find cliffs by enquing region information
  logmsg "Finding cliffs whose slopes are at least $$settings{slopeThreshold} with an rSquared value of at least $$settings{rSquaredThreshold} in a sliding window size of $$settings{windowSize}.";
  my $seqLength=bamSeqnames($bam,$settings);
  while(my($seqname,$length)=each(%$seqLength)){
    $Q->enqueue([$bam,$seqname,1,$length]);
  }

  # Join all the threads.
  $Q->enqueue(undef) for(@thr);
  my @cliff;
  for(@thr){
    my $cliff=$_->join();
    push(@cliff,@$cliff);
  }

  # Make regions contiguous if they are touching
  logmsg "Joining contiguous cliff regions";
  my %cliffRange;
  for my $cliff(@cliff){
    my $name=$$cliff{seqname};
    $cliffRange{$name}||=Array::IntSpan->new;
    # rounded to the nearest 5 or whatever the user specifies as a
    # minimum threshold.
    my $roundedSlope=int($$cliff{slope}/$$settings{slopeThreshold})*$$settings{slopeThreshold};
    $cliffRange{$name}->set_range($$cliff{start},$$cliff{stop},$roundedSlope);
  }

  # Print regions where cliffs are.
  # Include slope information as the score, but its directionality
  # as the strand.
  my $bedOut=Bio::FeatureIO->new(-format=>"bed"); # STDOUT
  for my $seqname(sort {$a cmp $b} keys(%cliffRange)){
    $cliffRange{$seqname}->consolidate(); # merge intersecting regions

    for my $range($cliffRange{$seqname}->get_range_list()){
      my($lo,$hi)=@$range;
      my $slope=$cliffRange{$seqname}->lookup($lo);
      my $score=abs($slope);

      # Positive slope: 1; Negative slope: 0; unknown: -1
      my $strand=1;
         $strand=0 if($slope < 0);

      my $name="$seqname:$lo";
      my $feat=Bio::SeqFeature::Annotated->new(-seq_id=>$seqname,-start=>$lo, -end=>$hi, -primary=>"Cliff", -Name=>$name, -source=>basename($0), -score=>$score, -strand=>$strand);
      $feat->annotation->add_Annotation("Name",Bio::Annotation::SimpleValue->new(-tagname=>"name",-value=>$name));
      $bedOut->write_feature($feat);
    }
  }

  return 0;
}

sub findCliffsBySlope{
  my($bam,$seqname,$pos,$lastPos,$settings)=@_;

  # Find the depth of the whole bam file one time, to save cpu.
  logmsg "Calculating depth of $bam:$seqname";
  my (undef,$depthfile)=tempfile("depthXXXXXX",SUFFIX=>".depth",DIR=>$$settings{tempdir});
  system("samtools depth -r '$seqname' $bam > $depthfile");
  die "ERROR with samtools depth: $!" if $?;
  my %depth;
  open(DEPTH,$depthfile) or die $!;
  while(<DEPTH>){
    chomp;
    my($seqname,$pos,$d)=split /\t/;
    print "$_\n" if(!$pos);
    $depth{$pos}=$d;
  }
  close DEPTH;


  my $windowSize=10;

  # Figure out the slope for every 15 bp in a sliding window
  # If it goes beyond the threshold then it is a cliff.
  my @cliff=();
  for($pos=$pos;$pos<=$lastPos;$pos++){

    # Load up information for line-fitting (ie, trendline)
    my (@depth,@pos); # y and x on the linefit graph
    my $start=$pos;
    my $stop=$start+$windowSize-1;
    for(my $windowPos=$start;$windowPos<=$stop;$windowPos++){
      my $depth=$depth{$windowPos}||0;
      push(@depth,$depth);
      push(@pos,$windowPos);
    }

    # Get the actual linefit statistics
    my $linefit=Statistics::LineFit->new();
    $linefit->setData(\@pos,\@depth) or die "Invalid depths to figure out the slope! \n @depth\n";

    # Don't worry about this being a cliff if rSquared is not significant
    next if(!defined($linefit->rSquared()) || $linefit->rSquared() < $$settings{rSquaredThreshold});

    # Get the slope of the line
    my($intercept,$slope)=$linefit->coefficients();

    # Don't think about this as a cliff if the slope is shallow
    next if(abs($slope) < $$settings{slopeThreshold});

    # Record this as a cliff if it passed all the criteria
    push(@cliff,{
      seqname=>$seqname, start=>$start, stop=>$stop,
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

# Use samtools view -H to find the seqnames and their lengths
sub bamSeqnames{
  my($bam,$settings)=@_;
  
  my %seqLength;

  open(BAMINFO,"samtools view -H '$bam' | ") or die "ERROR: could not use samtools view on $bam: $!";
  while(<BAMINFO>){
    chomp;
    my($type,@theRest)=split(/\t/,$_);
    next if($type ne '@SQ');
    my $theRest=join("\t",@theRest);

    my ($seqname,$length);
    while($theRest=~/(\w+?):([^\t]+)/g){
      $seqname=$2 if($1 eq "SN");
      $length=$2 if($1 eq "LN");
    }

    $seqLength{$seqname}=$length;
  }
  close BAMINFO;
  return \%seqLength;
}

sub usage{
  "Finds cliffs in the depth in a given sorted/indexed bam file.
  Usage: $0 file.bam
  --numcpus              1
  --tempdir              ''
  --windowSize           10   How many positions are used when calculating a slope
  --slopeThreshold       5    Ie, the depth changes X per every position in the window
  --rSquaredThreshold    0.7  The level of significance in the slope trendline
  "
}
