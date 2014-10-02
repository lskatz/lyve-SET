#!/usr/bin/env perl

# Changes a SNP matrix to an MSA
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use Bio::AlignIO;
use File::Basename qw/fileparse/;

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

  my $matrix=readMatrix($settings);
  $matrix=filterMatrix($matrix,$settings);

  printSeqs($matrix,$settings);

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
  my %matrix;
  while(<$fh>){
    chomp;
    my($genome,@nt)=split(/\t/,$_);
    for(my $i=0;$i<@nt;$i++){
      my ($contig,$pos)=split(/:/,$pos[$i]);
      $matrix{$genome}{$contig}{$pos}=$nt[$i];
    }
  }
  
  close $fh;
  return \%matrix;
}

sub filterMatrix{
  my($matrix,$settings)=@_;
  while(my($genome,$contigHash)=each(%$matrix)){
    while(my($contig,$posHash)=each(%$contigHash)){
      my $filteredHash=filterSortedPos($posHash,$settings);
      $$matrix{$genome}{$contig}=$filteredHash;
    }
  }
  return $matrix;
}

sub filterSortedPos{
  my($posHash,$settings)=@_;
  my @sorted=sort {$a<=>$b} keys(%$posHash);

  my %filtered;
  my $numPos=@sorted;

  for(my $i=1;$i<$numPos;$i++){
    my $pos=$sorted[$i];
    my $prevPos=$sorted[$i-1];
    my $distance=$pos - $prevPos;
    if($distance >= $$settings{'min-distance'}){
      $filtered{$pos}=$$posHash{$pos};
    } else {
      #logmsg "filtered: $pos is too close to $prevPos";
    }
  }
  return \%filtered;
}


sub printSeqs{
  my($matrix,$settings)=@_;


  my @seqObj;
  while(my($genome,$contigHash)=each(%$matrix)){
    my $sequence="";
    while(my($contig,$posHash)=each(%$contigHash)){
      # Create a list of sorted positions
      # The sorted positions assume that all sequences should have the same positions which is true! 
      # (but keep an eye on this just in case)
      my @sortedPos=sort {$a<=>$b} keys(%$posHash);
      for(@sortedPos){
        $sequence.=$$posHash{$_} || die "=>$_\n".Dumper [$contig,$posHash];
      }
    }
    my $seqObj=Bio::LocatableSeq->new(-id=>$genome,-seq=>$sequence,-start=>1,-end=>length($sequence));
    push(@seqObj,$seqObj);
  }

  my $out=Bio::AlignIO->new(-format=>$$settings{format},-displayname_flat => 1);
  my $aln=Bio::SimpleAlign->new(-seqs=>\@seqObj,-source=>"Lyve-SET:$0");
  $out->write_aln($aln);


  return 1;
}

sub usage{
  "$0: turns a SNP matrix into an alignment. The matrix can be the first parameter or stdin.
  Usage: $0 matrix.tsv > out.aln
         $0 < matrix.tsv | removeUninformativeSites.pl > informative.aln.fas
  --min-distance 0 The distance in bp between allowed SNPs. Does not count sites where all nt are the same.
  --format fasta   The output format (uses bioperl's aln format values, e.g., clustalw)
  "
}
