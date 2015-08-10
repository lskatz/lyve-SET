#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help));
  $$settings{maskedThresholdPercent}=10;

  die usage if($$settings{help} || !@ARGV);

  for my $matrix(@ARGV){
    reportMaskedGenomes($matrix,$settings);
  }
}

sub reportMaskedGenomes{
  my($matrix,$settings)=@_;
  open(MATRIX,$matrix) or die "ERROR: could not open $matrix for reading: $!";
  while(<MATRIX>){
    

  close MATRIX;
}
