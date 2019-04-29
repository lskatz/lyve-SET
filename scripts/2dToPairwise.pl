#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help));
  die usage() if($$settings{help});

  # read the file
  my $header = <STDIN>;
  if(!$header){
    print STDERR "$0: ERROR: not stdin was given\n";
    return 1;
  }

  chomp($header);
  my @header = split(/\t/, $header);
  my $numFields = @header;
  my %matrix;
  while(<STDIN>){
    s/^\s+|\s+$//g; # trim whitespace
    next if(/^$/);  # skip empty lines
    my @F=split /\t/;
    die "ERROR: this line does not have $numFields fields\n".usage() if(@F<$numFields);

    my $ref = $F[0];
    for(my $i=1;$i<$numFields; $i++){
      my $value = $F[$i];
      my $query = $header[$i];

      my($g1, $g2) = sort($ref, $query);

      $value = 0 if($value eq '-');
      $matrix{$g1}{$g2} = $value;
    }
  }
  close STDIN;

  while(my($g1, $distances) = each(%matrix)){
    while(my($g2, $dist) = each(%$distances)){
      print join("\t", $g1, $g2, $matrix{$g1}{$g2})."\n";
    }
  }
  return 0;
}

sub usage{
  $0=basename $0;
  "Turns a 2d matrix/spreadsheet into a 3-column format
  Usage: $0 < table.tsv > pairwise.tsv
  "
}
