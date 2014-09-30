#!/usr/bin/env perl

# Changes a SNP matrix to a MSA
# Author: Lee Katz

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use File::Basename;

$0=fileparse $0;
sub logmsg{$|++;print STDERR "$0: @_\n"; $|--;}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help verbose min-distance=i format=s)) or die $!;
  die usage() if($$settings{help});
  $$settings{"min-distance"} ||=0;
  $$settings{format}         ||="fasta";

  die usage() if($$settings{help});

  my $seqObj=readMatrix($settings);

  printSeqs($seqObj,$settings);

  return 0;
}

sub readMatrix{
  my ($settings)=@_;

  # read either stdin or the filename
  my $fh;
  if(@ARGV == 1){
    logmsg "Reading from $ARGV[0]";
    open($fh,$ARGV[0]);
  } elsif(@ARGV > 1){
    die "ERROR: you can only have one filename but you gave multiples: ".join(", ",@ARGV)."\n".usage();
  } else {
    $fh=\*STDIN;
    logmsg "Reading from stdin";
  }

  my $posHeader=<$fh>;
  chomp($posHeader);
  my ($dot,@pos)=split(/\t/,$posHeader);
  my $seqObj;
  while(<$fh>){
    chomp;
    my($genome,@nt)=split(/\t/,$_);
    my $sequence=join("",@nt);
    my $seq=Bio::Seq->new(-id=>$genome,seq=>$sequence);
    push(@$seqObj,$seq);
  }
  
  close $fh;
  return ($seqObj,\@pos) if wantarray;
  return $seqObj;
}

sub printSeqs{
  my($seq,$settings)=@_;
  my $out=Bio::SeqIO->new(-format=>$$settings{format},-fh=>\*STDOUT);
  return $out->write_seq(@$seq);
}

sub usage{
  "$0: turns a SNP matrix into an alignment
  Usage: $0 < matrix.tsv > aln.out.fas
         $0 < matrix.tsv | removeUninformativeSites.pl > informative.aln.fas
  --min-distance 0 The distance in bp between allowed SNPs. Does not count sites where all nt are the same.
  --format fasta   The output format
  "
}
