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

      $value = 0 if($value eq '-');
      $matrix{$ref}{$query} = $value;
    }
  }
  close STDIN;

  my @sortedRef = sort{$a cmp $b} grep {$_ ne '.'} values(@header);
  for(my $i=0; $i<@sortedRef; $i++){
    for(my $j=0; $j<@sortedRef; $j++){
      my $query = $sortedRef[$j] || "";
      my $ref = $sortedRef[$i] || "";
      print join("\t", $ref,
                       $query,
                       $matrix{$ref}{$query},
                     );
      print "\n";
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
