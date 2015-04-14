#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use List::Util qw/sum/;
use Statistics::Descriptive;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;

exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help));
  
  my $bam=shift(@ARGV);
  die "ERROR: bam file not given!\n".usage() if(!$bam);
  die "ERROR: could not find $bam\n".usage() if(!-e $bam);

  my @region=findCliffs($bam,$settings);

  return 0;
}

sub findCliffs{
  my($bam,$settings)=@_;
  my %region; # $region{contig}=>[ 
              #                    [start,stop],
              #                    [start,stop],
              #                  ]

  my @seqname=getSeqname($bam,$settings);
  my $windowSize=20;
  my $stepSize=10;
  my $numDataPoints=100;
  my $lastPos=$windowSize*$numDataPoints;
  for my $seqname(@seqname){
    # Find the average of average depths
    my %depth;
    for(my $start=1;$start<$lastPos;$start+=$stepSize){
      my $stop=$start+$windowSize-1;
      $stop=$lastPos if($stop>$lastPos);
      my $windowDepth=depth($bam,"$seqname:$start-$stop",$settings);
      $depth{$start}=$windowDepth;
    }
    my $stat=Statistics::Descriptive::Full->new;
    $stat->add_data(values(%depth));

    # Is anything lower than the 5% cutoff?
    my $lo=$stat->mean - 2*$stat->standard_deviation;
    my $hi=$stat->mean + 2*$stat->standard_deviation;
    while(my($start,$depth)=each(%depth)){
      next if($depth<$hi && $depth>$lo);
      push(@{$region{$seqname}},[$start,$start+$windowSize]);
    }
    #logmsg join("\t",$stat->mean,$stat->standard_deviation,$lo,$hi,Dumper(\%region))."\n";
  }
  die;
}

sub getSeqname{
  my($bam,$settings)=@_;
  my (@seqname);
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
    push(@seqname,$F{SN}) if($F{SN});
  }
  return @seqname;
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
  
  return sum(@depth)/scalar(@depth);
}

sub usage{

}
