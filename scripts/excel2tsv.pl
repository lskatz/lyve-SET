#!/usr/bin/env perl

use strict;
use warnings;
no warnings 'utf8';  #hides 'Wide character in print' error
use File::Basename qw/basename fileparse/;
use Getopt::Long;

use FindBin;
use lib "$FindBin::RealBin/../lib/lib/perl5";
use Spreadsheet::ParseExcel;
use Spreadsheet::XLSX;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;

exit(main());

sub main{

  # Get any command-line settings
  my $settings={};
  GetOptions($settings,qw(prefix=s help));
  die usage() if($$settings{help});
  my $file=$ARGV[0];
  die "ERROR: need Excel file\n".usage() if(!$file);
  $$settings{prefix}||="$file";

  # Figure out path and extension information
  my($name,$path,$extension)=fileparse($file,qw(.xlsx .xls));
  $extension=lc($extension); # lowercase to make comparisons easier

  if($extension eq '.xls'){

      my $parser = Spreadsheet::ParseExcel -> new();
      my $workbook = $parser -> parse($file);
      if(!defined $workbook){
          die $parser -> error(), ".\n";
      }

      for my $sheet($workbook -> worksheets()){
          my $sheetName = $sheet -> get_name();
          my $filename = "$$settings{prefix}"."_$sheetName.tsv";
          my($row_min,$row_max)=$sheet -> row_range();
          my($col_min,$col_max)=$sheet -> col_range();
          open(OUTFILE,'>',$filename) or die "ERROR: could not open $filename for writing: $!";
          for my $row ($row_min .. $row_max){
              my @values;
              for my $col ($col_min .. $col_max){
                  my $cell = $sheet -> get_cell($row, $col);
                  if(defined $cell and $cell -> value() ne ''){
                      push @values, $cell -> value();
                  }
                  else{
                      push @values, "\t";
                  }
              }
              $" = "\t";
              print OUTFILE "@values\n";
          }
          close OUTFILE;
      }
      logmsg "XLS converted into TSV\n";
  }

  elsif($extension eq '.xlsx'){

      my $workbook = Spreadsheet::XLSX -> new ($file);
      foreach my $sheet (@{$workbook -> {Worksheet}}){
          my $sheetName = $sheet->{Name};
          my $filename = "$$settings{prefix}"."_$sheetName.tsv";
          open(OUTFILE,'>',$filename) or die "ERROR: could not open $filename for writing: $!";
          $sheet -> {MaxRow} ||= $sheet -> {MinRow};    
          foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow}){             
              $sheet -> {MaxCol} ||= $sheet -> {MinCol};                
              foreach my $col ($sheet -> {MinCol} .. $sheet -> {MaxCol}){
                  my $cell = $sheet -> {Cells} [$row] [$col];                
                  if($cell){
                      print OUTFILE $cell -> {Val}, "\t";
                  }
                  else{
                      print OUTFILE "\t\t";
                  }
              }
              print OUTFILE "\n";
          }
      }
      logmsg "XLSX converted into TSV\n";
  }

  else{
      print "ERROR: unsupported extension\n";
      print "This script requires a XLS or XLSX extension but you supplied $extension\n\n";
      usage();
  }

  return 0;
}

sub usage{
local $0=basename($0);
"Converts a Microsoft Excel spreadsheet into TSV format.
Output is automatically created in same directory as the input.
Each worksheet is written as a separate file.

Usage: $0 <inputfile>
  -p outputPrefix
";
}
