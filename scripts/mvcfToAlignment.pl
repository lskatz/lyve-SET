#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use File::Basename;

$0=fileparse $0;

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help));
  die usage() if($$settings{help});
  
  my $header=<>;
  chomp($header);
  $header=~s/^#\s*//; # remove the hash in front of the header
  my @header=split(/\t/,$header);
  $_=~s/\[\d+\]// for(@header); # remove [number] notations for the headers

  # genome names
  my @genome=@header[3..@header-1];

  # read the matrix into memory before printing the aln
  my %matrix;
  while(<>){
    chomp;
    my %F;
    @F{@header}=split /\t/;

    # capture and remove non-nucleotide information from the matrix
    my %extra;
    for(qw(CHROM  POS ID  REF  ALT QUAL FILTER  INFO    FORMAT)){
      $extra{$_}=$F{$_};
      delete($F{$_});
    }

    # assign the nts for each genome, at each contig/pos
    while(my($genome,$GT)=each(%F)){
      die Dumper [$genome,\%F] if(!$GT);
      my $nt=$GT;
      if($nt=~/(.)[\/\|](.)/){
        $nt=$1;
        $nt="N" if($1 ne $2);
        $nt=$extra{REF} if($nt eq '.' || $nt eq ',');
        
        # take care of indels or ambiguities
        $nt="N" if(length($nt)!=1 || $nt!~/^[ATCG]$/i);
      }
      $matrix{$genome}{$extra{CHROM}}{$extra{POS}}=$nt;
    }
  }
  # DONE getting it into memory!

  # START output
  # get all the sorted positions
  my @CHROM=keys(%{ $matrix{$genome[0]} });
  my @POS;
  for my $chr(sort{$a cmp $b} @CHROM){
    push(@POS,"$chr:$_") for(sort{$a<=>$b} keys(%{ $matrix{$genome[0]}{$chr} }));
  }

  # create the alignment
  my $out=Bio::SeqIO->new(-format=>"fasta");
  for my $genome(@genome){
    my $seq="";
    for(@POS){
      my($chr,$pos)=split(/:/,$_);
      $seq.=$matrix{$genome}{$chr}{$pos};
    }
    my $seqObj=Bio::Seq->new(-seq=>$seq,-id=>$genome);
    $out->write_seq($seqObj);
  }
  
  return 0;
}

sub usage{
  " $0: reads a bcftools query output and creates a multiple sequence alignment file
  Usage: bcftools query [...] | $0 > aln.fasta
  "
}
