#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Fixes and adds fields on a given Vcf

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename fileparse/;
use File::Temp qw/tempdir/;
use Bio::Perl;
use Bio::FeatureIO;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/@fastaExt @fastqExt @bamExt/;
#use Vcf;

$ENV{PATH}="$FindBin::RealBin:$ENV{PATH}";
$ENV{BCFTOOLS_PLUGINS}="$FindBin::RealBin/../lib/bcftools-1.2/plugins";
$ENV{LD_LIBRARY_PATH}=$ENV{LD_LIBRARY_PATH}.":$FindBin::RealBin/../lib/bcftools-1.2/htslib-1.2.1";
use constant reportEvery => 1000000;

$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());
 

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s numcpus=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $VCF=$ARGV[0];
  die usage() if(!$VCF || $$settings{help});
  
  my $fixedVcf=fixVcf($VCF,$settings);

  logmsg "Done fixing it up! Printing to stdout.";
  system("cat $fixedVcf");
  die if $?;

  return 0;
}

sub fixVcf{
  my($vcf,$settings)=@_;
  
  my $uncompressed=bcftoolsFixes($vcf,$settings);
  my $richVcf=addFeatureTags($uncompressed,$settings);

  return $richVcf;
}

sub bcftoolsFixes{
  my($compressedVcf,$settings)=@_;
  my $v="$$settings{tempdir}/filledInAnAcRef.vcf";
  logmsg "Filling in AN and AC, and also missing REF in the vcf";
  
  #system("LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$FindBin::RealBin/../lib/bcftools-1.2/htslib-1.2.1 bcftools plugin fill-AN-AC $compressedVcf | bcftools plugin missing2ref > $v");
  system("bcftools plugin fill-AN-AC $compressedVcf | bcftools plugin missing2ref > $v");
  die "ERROR: could not run bcftools plugins on $compressedVcf. This is what plugins are available:\n `bcftools plugin -lv 2>&1`:\n".`bcftools plugin -lv 2>&1` if $?;

  return $v;
}

# add FT tag
sub addFeatureTags{
  my($vcf,$settings)=@_;
  my $v="$$settings{tempdir}/addedFt.vcf";

  logmsg "Adding feature tags that are missing";

  my $numSites=0;
  my $hasFtTag=0; # Whether the input file already has the Ft tag in the header or whether it has already been added to the output
  open(VCFIN,$vcf) or die "ERROR: could not read $vcf: $!";
  open(VCFOUT,">", $v) or die "ERROR: could not write to $v: $!";
  while(<VCFIN>){
    if(/^##/){
      $hasFtTag=1 if(/FORMAT=<ID=FT/);
      print VCFOUT $_;
      next;
    }
    elsif(/^#CHROM/){
      # Print extra headers and then the column labels last
      print VCFOUT '##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filters using the same codes as the FILTER data element">'."\n" if(!$hasFtTag);
      print VCFOUT $_;
      next;
    }

    if(++$numSites % reportEvery == 0){
      logmsg "Looked at $numSites sites";
    }

    # Read the Vcf line
    my($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@SAMPLE)=split(/\t/);
    chomp(@SAMPLE);

    # Does it already have the FT tag? If not, then add it
    my(@FORMAT)=split(/:/,$FORMAT);
    my $hasFtValue=0;
       $hasFtValue=1 if(grep(/^FT$/,@FORMAT));
    if(!$hasFtValue){
      push(@FORMAT,'FT');
    }
    $FORMAT=join(':',@FORMAT);

    # Look at each sample to see if the FT tag has a value
    # or if other tags have values.
    my $numSamples=@SAMPLE;
    for(my $i=0; $i<$numSamples; $i++){
      my %sampleInfo;
      @sampleInfo{@FORMAT}=split(/:/,$SAMPLE[$i]);

      # Fill in any specific undefined values with what I want
      $sampleInfo{'FT'}=$FILTER if(!defined($sampleInfo{FT}));

      # Fill in the rest of the undefined values
      while(my($key,$value)=each(%sampleInfo)){
        $sampleInfo{$key}='.' if(!defined($sampleInfo{$key}));
      }

      # join the line back together in the proper order
      $SAMPLE[$i]='';
      for my $f(@FORMAT){
        $SAMPLE[$i].="$sampleInfo{$f}:";
      }
      $SAMPLE[$i]=~s/(:\.)*:$//; # remove last semicolon and also any blank values just hanging out at the end of the line
    }

    my $modifiedLine=join("\t",$CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@SAMPLE);
    print VCFOUT "$modifiedLine\n";
  }
  close VCFOUT;
  close VCFIN;

  return $v;
}

sub usage{
  "Fixes a given VCF by adding, removing, or editing fields
  Usage: $0 file.vcf.gz > fixed.vcf
  --numcpus 1
  "
}
