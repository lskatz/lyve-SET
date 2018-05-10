#!/usr/bin/env perl

require 5.12.0;
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use File::Basename;
use File::Temp qw/tempdir/;
use List::Util qw/min max sum uniq/;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;
use lib "$FindBin::RealBin/../lib/lib/perl5";
use Number::Range;

use constant reportEvery=>100000;

$0=fileparse $0;

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s density-window|window=i density-filter=i)) or die $!;
  die usage() if($$settings{help});
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1,CLEANUP=>1);

  # Ensure that the density filter and density window are >= 1
  if(!$$settings{'density-filter'} || $$settings{'density-filter'} < 1){
    $$settings{'density-filter'}=1;
  }
  if(!$$settings{'density-window'} || $$settings{'density-window'} < 1){
    $$settings{'density-window'}=1;
  }

  my($in)=@ARGV;
  $in or die "ERROR: need input file\n".usage();
  if(!-e $in){
    die "ERROR: could not find input file $in";
  }

  my $sites = sitesToInclude($in,$settings);

  if(-s $sites == 0){
    die "ERROR: all sites are excluded using the density filter of ".$$settings{'density-filter'}." in a window size ".$$settings{'density-window'};
  }

  system("bcftools filter -R $sites '$in'"); die if $?;

  return 0;
}

sub sitesToInclude{
  my($vcf,$settings)=@_;

  # Open the vcf
  open(my $fp,"zcat $vcf | ") or die "ERROR: could not zcat $vcf: $!";
  <$fp>; # skip the header

  # save all positions in a hash of arrays
  my %pos;
  while(my $vcfLine=<$fp>){
    next if($vcfLine=~/^#/);
    # Save the positional information
    my($seq,$pos)=split(/\t/,$vcfLine);
    push(@{$pos{$seq}}, $pos);
  }
  close $fp;

  # run filtering on each array and mark positions
  my $includedSites="$$settings{tempdir}/regions.tsv";
  open(my $includedSitesFh, ">", $includedSites) or die "ERROR writing to $includedSites: $!";
  while(my($seq,$posArr) = each(%pos)){
    my $filteredPos = densityFilter($posArr,$$settings{'density-filter'},$$settings{'density-window'},$settings);
    for my $pos(@$filteredPos){
      print $includedSitesFh join("\t",$seq,$pos)."\n";
    }
  }
  close $includedSitesFh;
  
  return $includedSites;
}

# For an array of positions, run the density filter
sub densityFilter{
  my($posArr,$densityFilter,$densityWindow, $settings)=@_;

  # We will return filtered sites
  my @filtered;

  my @pos = sort {$a <=> $b} @$posArr;
  my $numPos=@pos;
  my $windowStartIndex=0;
  my $lastPos=$numPos - $densityFilter; # To capture the max value in the for loop
  for(my $i=0;$i<$lastPos;$i++){
    # Advance ahead the number of sites allowed in the window.
    # Is the window size greater than the density-filter-window?
    my $maxBound=$pos[$i+$densityFilter - 1];
    my $minBound=$pos[$i];
    my $window=$maxBound - $minBound;
    if($window >= $densityWindow){
      push(@filtered,$pos[$i]);
    }
  }
  
  # Capture the last X sites if they are not within the window size
  my $maxBound=$pos[$numPos-1];
  my $minBound=$pos[$numPos - $densityFilter - 1];
  my $window = $maxBound - $minBound;
  if($window >= $densityWindow){
    for(my $i=$numPos - $densityFilter;$i<$numPos;$i++){
      push(@filtered,$pos[$i]);
    }
  }

  # Sort/unique all positions just in case the two loops
  # are overlapping.
  @filtered = uniq(@filtered);
  @filtered = sort {$a <=> $b} @filtered;

  return \@filtered;
      
}

sub usage{
  "$0: filters a vcf file so that no window has too many sites. 
  The first three columns of the matrix are contig/pos/ref, and the next columns are all GT.

  A site will get filtered if it is upstream a set of sites within
  a given window size.  That window size is inclusive of start/stop
  positions.  For example if the window size is 10 and the filter
  density is 2, then there cannot be more than 2 sites between 
  positions 1 and 10. This includes positions 1 and 10.

  Usage: 
    $0 file.vcf.gz > filtered.vcf
    Default is to not filter.
  --density-filter    1         Number of sites allowed per window
  --density-window    1         Size of the window
  "
}

