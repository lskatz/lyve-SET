#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help outfile=s indels! depth=i badsites=s));
  my @infile=@ARGV;
  $$settings{outfile}||="$0.out.vcf";
  $$settings{depth}||=0;
  $$settings{indels}=1 if(!defined($$settings{indels}));
  die usage($settings) if($$settings{help} || !@infile);

  filterVcfs(\@infile,$settings);

  print "VCF output is in $$settings{outfile}\n";
  return 0;
}

sub filterVcfs{
  my($infile,$settings)=@_;

  my $newVcfStr="";
  open(OUT,">",$$settings{outfile}) or die "Could not open $$settings{outfile} for writing:$!";
  $newVcfStr.=filterVcf($_,$settings) for(@$infile);
  print OUT $newVcfStr;
  close OUT;
}

sub filterVcf{
  my($infile,$settings)=@_;
  my $vcfStr="";
  my @badSites;
  open(IN,"<",$infile) or die "Could not open $infile for reading:$!";
  while(<IN>){
    # print comments and then move on
    if(/^#/){
      $vcfStr.=$_;
      next;
    }

    # get some field information
    my @F=split /\t/;
    my $numFields=@F;
    my ($rseq,$pos)=@F[0..1];
    my $posId=join("_",$rseq,$pos);

    # find out if it's an indel
    if(!$$settings{indels} && (length($F[3]) !=  length($F[4])) ){
      push(@badSites,$posId);
      next;
    }

    # grab other info about the snp call
    my @info=split(/;/,$F[7]);
    my %info;
    for(@info){
      my($key,$value)=split /=/;
      $key=uc($key); # just to be sure
      $info{$key}=$value;
    }

    # If the cigar line has any indels indicated and indels are not wanted, move on.
    # CIGAR is defined in the SAM specifications, http://samtools.sourceforge.net/SAMv1.pdf
    $info{CIGAR}||=""; # make sure CIGAR is defined before testing it
    next if(!$$settings{indels} && $info{CIGAR}=~/[ID]/);

    # multibase snps are printed to a line per position, but this works for single base calls too
    my @alt=split(//,$F[4]);
    my @ref=split(//,$F[3]);
    my $numBases=max(scalar(@alt),scalar(@ref));
    for(my $i=0;$i<$numBases;$i++){
      $posId=join("_",$rseq,($pos+$i));
      # find out if it passes other attributes like depth
      if($info{DP}<$$settings{depth}){
        push(@badSites,$posId);
        next;
      }

      chomp(@F);
      $ref[$i]||="*";
      $alt[$i]||="*";
      next if($alt[$i] eq $ref[$i]); # if equal, then not a variant
      my $newLine=join("\t",$rseq,($pos+$i),$F[2],$ref[$i],$alt[$i],@F[5..($numFields-1)]);
      $vcfStr.="$newLine\n";
    }
  }
  close IN;
  reportBadSites(\@badSites,$settings) if($$settings{badsites});
  return $vcfStr;
}

sub reportBadSites{
  my($badSites,$settings)=@_;
  open(BADSITES,">",$$settings{badsites}) or die "Could not open $$settings{badsites} for writing: $!";
  print BADSITES $_."\n" for(@$badSites);
  close BADSITES;
  print "Bad sites have been printed to $$settings{badsites}\n";
  return 1;
}

sub usage{
  "Filters out variant calls that you do not want in a VCF file
  usage: $0 infile.vcf [-o outfile.vcf]
  --noindels to remove indels
  -d depth
    to remove variants with this depth and lower
  -b badsites.txt
    a file to remember where positions with low quality variations are held
  "
}
