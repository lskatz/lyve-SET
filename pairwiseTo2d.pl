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

  # spit out the spreadsheet
  my @header=sort {$a cmp $b} keys(%matrix);
  my $numHeaders=@header;
  print join("\t",".",@header)."\n";
  for(my $i=0;$i<$numHeaders;$i++){
    my $key1=$header[$i];
    print "$key1\t";
    for(my $j=0;$j<$numHeaders;$j++){
      my $key2=$header[$j];
      my $value=$matrix{$key1}{$key2} || '-';
      print $value."\t";
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
