#!/usr/bin/env perl

require 5.12.0;
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use File::Basename;
use File::Temp qw/tempdir/;
use List::Util qw/min max sum/;

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

  $$settings{'density-filter'}||=0;
  $$settings{'density-window'}||=0;

  my($in)=@ARGV;
  $in or die "ERROR: need input file\n".usage();

  my $sites = sitesToInclude($in,$settings);
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
    my $filteredPos = densityFilter($posArr,$settings);
    for my $pos(@$filteredPos){
      print $includedSitesFh join("\t",$seq,$pos)."\n";
    }
  }
  close $includedSitesFh;
  
  return $includedSites;
}

# For an array of positions, run the density filter
sub densityFilter{
  my($posArr, $settings)=@_;

  my @filtered;

  my @pos = sort {$a <=> $b} @$posArr;
  my $numPos=@pos - $$settings{'density-filter'};
  my $windowStartIndex=0;
  for(my $i=0;$i<$numPos;$i++){
    # Advance ahead the number of sites allowed in the window.
    # Is the window size greater than the density-filter-window?
    my $window=$pos[$i+$$settings{'density-filter'}] - $pos[$i];
    if($window > $$settings{'density-window'}){
      push(@filtered,$pos[$i]);
    }
  }
  
  # Capture the last X sites if they are not within the window size
  my $window = $pos[$numPos + $$settings{'density-filter'} - 1] - $pos[$numPos];
  if($window > $$settings{'density-window'}){
    for(my $i=$numPos;$i<@pos;$i++){
      push(@filtered,$pos[$i]);
    }
  }

  return \@filtered;
      
}

sub usage{
  "$0: filters a vcf file so that no window has too many sites. 
  The first three columns of the matrix are contig/pos/ref, and the next columns are all GT.

  Usage: 
    $0 file.vcf.gz > filtered.vcf
    Default is to not filter.
  --density-filter    1         Number of sites allowed per window
  --density-window    1         Size of the window
  "
}
