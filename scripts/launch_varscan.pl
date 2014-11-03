#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Creates a vcf using varscan

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Bio::Perl;

$0=fileparse $0;
sub logmsg{print STDERR "@_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help reference=s tempdir=s altFreq=s coverage=i));

  my $bam=$ARGV[0] or die "ERROR: need bam\n".usage();
  my $reference=$$settings{reference} or die "ERROR: need --reference\n".usage();
  $$settings{tempdir}||="tmp";
  $$settings{coverage}||=10;
  $$settings{altFreq}||=0.75;

  # Check to see if varscan is installed correctly
  `varscan.sh >/dev/null`;
  die "ERROR: varscan.sh gave an error" if $?;

  my $pileup=mpileup($bam,$reference,$settings);
  varscan($pileup,$settings);
  return 0;
}

sub mpileup{
  my($bam,$reference,$settings)=@_;
  my $pileup="$$settings{tempdir}/".fileparse($bam).".mpileup";
  return $pileup if(-e $pileup && -s $pileup > 0);
  logmsg "Creating a pileup $pileup";
  system("samtools mpileup -f '$reference' '$bam' 1>$pileup");
  die if $?;
  return $pileup;
}

sub varscan{
  my($pileup,$settings)=@_;
  die "ERROR: the pileup is a zero-byte file\n  $pileup" if(-s $pileup < 1);
  system("varscan.sh mpileup2cns $pileup --min-coverage $$settings{coverage} --min-coverage 10 --min-var-freq 0.75 --output-vcf 1");
  die if $?;
}

sub usage{
  local $0=fileparse $0;
  "$0: find SNPs using varscan
  Usage: $0 file.bam > file.vcf
  -t tmp/ A temporary directory to store files
  -c 10   Min coverage
  -a 0.75 Min consensus agreement for a SNP
  "
}
