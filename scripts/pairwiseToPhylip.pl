#!/usr/bin/env perl

# turns a 3-column pairwise distance into a 2d matrix/spreadsheet

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help));
  die usage() if($$settings{help});

  # read the file
  my %matrix;
  while(<>){
    s/^\s+|\s+$//g; # trim whitespace
    next if(/^$/);  # skip empty lines
    my @F=split /\t/;
    die "ERROR: this line does not have three fields\n".usage() if(@F<3);
    $matrix{$F[0]}{$F[1]}=$F[2];
    $matrix{$F[1]}{$F[0]}=$F[2];
  }

  my @taxon   = keys(%matrix);
  my $numTaxa = @taxon;

  print "    $numTaxa\n";
  for(my $i=0;$i<$numTaxa;$i++){
    my $refTaxon = $taxon[$i];
    print $refTaxon;
    for(my $j=0;$j<$numTaxa;$j++){
      my $queryTaxon = $taxon[$j];
      my $dist = $matrix{$refTaxon}{$queryTaxon};
      if(!defined($dist)){
        if($i == $j){
          $dist = 0;
        } else {
          die "ERROR: I cannot find the distance between $refTaxon and $queryTaxon";
        }
      }
      print "\t" . $dist;
    }
    print "\n";
  }
    
  return 0;
}

sub usage{
  $0=fileparse $0;
  "Turns a 3-column pairwise distance into a 2d matrix/spreadsheet
  Usage: $0 < pairwise.tsv > table.tsv
  "
}
