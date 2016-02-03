#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use Bio::AlignIO;
use Data::Dumper;
use Getopt::Long;
use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;

exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i header));
  die usage() if($$settings{help});
  $$settings{numcpus}||=1;

  die usage() if(@ARGV<1);

  for my $matrix(@ARGV){
    logmsg "Reading $matrix for pairwise SNPs";
    matrixDistances($matrix,$settings);
  }
  
  return 0;
}

sub matrixDistances{
  my($matrix,$settings)=@_;

  my %dist;
  my %denominator;

  open(MATRIX,"<",$matrix) or die "ERROR: could not open $matrix for reading: $!";
  my (@header,@posHeader,@sample,%header,$numSamples);
  while(<MATRIX>){
    chomp;
    
    # Grab the header
    if(/^#.*CHROM/){
      # clear the header
      $_=[] for(\@header,\@posHeader,\@sample);

      # Figure out the header values
      s/^#\s*//;
      @header=split /\t/;
      for(@header){
        $_=~s/\[\d+\]//;
        $_=~s/:GT$//;
      }
      @header{@header}=split(//,1 x scalar(@header)); # index the header

      # Grab potential columns that are not GT
      for my $h(qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT)){
        if($header{$h}){
          push(@posHeader,$h);
          delete($header{$h});
        }
      }
      @sample=keys(%header);
      $numSamples=scalar(@sample);

      next; # done getting headers!
    }

    # Get information on this line
    my %pos;
    @pos{@header}=split /\t/;

    # Make all the nts lowercase for ease of comparison
    for my $sample(@sample){
      $pos{$sample}=lc($pos{$sample});
    }

    for(my $i=0;$i<$numSamples;$i++){
      next if($pos{$sample[$i]} =~ /[^atcg]/);
      for(my $j=$i+1;$j<$numSamples;$j++){
        # Don't count ambiguities
        next if($pos{$sample[$j]} =~ /[^atcg]/);
      
        # Alphabetize just to make sure the hash is stable
        my($sample1,$sample2)=sort {$a cmp $b} ($sample[$i],$sample[$j]);

        # Now that this is a comparison, add onto the denominator.
        $denominator{$sample1}{$sample2}++;
        
        # Don't count same-nts as distance
        next if($pos{$sample1} eq $pos{$sample2});

        # Now that they are different and unambiguous, add the distance
        $dist{$sample1}{$sample2}++;
      }
    }
  }
  close MATRIX;

  # Print the distances
  if($$settings{header}){
    print join("\t",qw(GENOME1 GENOME2 SNP SITES FREQ))."\n";
  }
  for(my $i=0;$i<$numSamples;$i++){
    for(my $j=$i+1;$j<$numSamples;$j++){
      my ($sample1,$sample2)=sort {$a cmp $b} ($sample[$i],$sample[$j]);
      my $dist=$dist{$sample1}{$sample2} || 0;
      my $sites=$denominator{$sample1}{$sample2};
      if($sites < 1){
        logmsg "ERROR: no sites were comparable between $sample1 and $sample2. This could be caused by large numbers of ambiguous SNP calls.";
        next;
      }

      my $freq=sprintf("%0.2f",$dist/$sites);
      print join("\t",$sample1,$sample2,$dist,$sites,$freq)."\n";
    }
  }
}

sub usage{
  "Finds pairwise distances of entries in a SNP matrix.
  SNP matrix format must be at least three fields: CHROM, REF, ALT, [GT, ...]
  and has an optional header with a leading hash sign (#).
  Usage: $0 snpmatrix.tsv > pairwise.tsv
  --header          Display a header for the columns        
  --numcpus   1
  "
}
