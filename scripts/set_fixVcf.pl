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
use Vcf;

$ENV{PATH}="$FindBin::RealBin:$ENV{PATH}";
$ENV{BCFTOOLS_PLUGINS}="$FindBin::RealBin/../lib/bcftools-1.2/plugins";
$ENV{LD_LIBRARY_PATH}||="";
$ENV{LD_LIBRARY_PATH}=$ENV{LD_LIBRARY_PATH}.":$FindBin::RealBin/../lib/bcftools-1.2/htslib-1.2.1";
use constant reportEvery => 1000000;

$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());
 

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s numcpus=i min_coverage=i min_alt_frac=s fail-samples fail-sites rename-sample=s filter-ft add-N removed-unused-alt)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1,CLEANUP=>1);
  mkdir $$settings{tempdir};
  $$settings{min_coverage}||=0;
  $$settings{min_alt_frac}||=0;
  die "ERROR: min_alt_frac must be between 0 and 1, inclusive" if($$settings{min_alt_frac} < 0 || $$settings{min_alt_frac} > 1);

  my $VCF=$ARGV[0];
  die usage() if(!$VCF || $$settings{help});
  
  my $fixedVcf=fixVcf($VCF,$settings);

  logmsg "Printing to stdout.";
  system("cat $fixedVcf");
  die if $?;

  return 0;
}

sub fixVcf{
  my($vcf,$settings)=@_;
  
  logmsg "Fixing the VCF";
  my $uncompressed=addStandardTags($vcf,$settings);
  #my $richVcf=addFilterTag($uncompressed,$settings);
  my $reevaluated=reevaluateSites($uncompressed,$settings);
  logmsg "Done fixing it up!";

  return $reevaluated;
}

sub addStandardTags{
  my($compressedVcf,$settings)=@_;
  my $v="$$settings{tempdir}/filledInAnAcRef.vcf";
  logmsg "Filling in AN and AC, and also missing REF in the vcf";
  
  system(qq(
    bcftools view $compressedVcf --exclude-types indels |\
    bcftools plugin fill-AN-AC > $v
    )
  );
  die "ERROR: could not run bcftools plugins on $compressedVcf. This is what plugins are available:\n `bcftools plugin -lv 2>&1`:\n".`bcftools plugin -lv 2>&1` if $?;

  # rename the sample if requested
  if(my $newname=$$settings{'rename-sample'}){
    my $samplesFile="$$settings{tempdir}/samples.txt";

    # Bcftools reheader is weird. Need to bgzip and need to
    # account for bgzip'd output.
    # Also need to make a samplename file.
    system("echo '$newname' > $samplesFile"); die if $?;
    system("bgzip $v && tabix $v.gz");
    system("bcftools reheader --samples $samplesFile $v.gz | zcat > $v");
    die "ERROR with renaming the sample" if $?;
    system("rm -f $v.gz $v.tbi");  # cleanup
  }

  return $v;
}

sub reevaluateSites{
  my($vcf,$settings)=@_;
  my $v="$$settings{tempdir}/reevaluated.vcf";
  logmsg "Reevaluating sites on whether they pass thresholds";

  my $altFreq=$$settings{min_alt_frac};
  my $inverseAltFreq=1-$altFreq;

  my $vcfObj=Vcf->new(file=>$vcf);
  open(OUTVCF,">",$v) or die "ERROR: could not open vcf file $v for writing: $!";
  
  $vcfObj->parse_header(); # I might need to parse headers, first thing...?
  my (@samples) = $vcfObj->get_samples(); # have to parse headers before getting samples
  my $numSamples=scalar(@samples);
  $vcfObj->add_header_line({key=>'FILTER', ID=>"FAIL", Description=>"This site did not pass QC"}) if(!@{ $vcfObj->get_header_line(key=>'FILTER',ID=>'FAIL') });
  $vcfObj->add_header_line({key=>'FORMAT', ID=>'FT', Number=>1, Type=>'String', Description=>"Genotype filters using the same codes as the FILTER data element"}) if(!@{ $vcfObj->get_header_line(key=>'FORMAT',ID=>'FT') });

  # start printing
  my $vcfHeader=$vcfObj->format_header();
  print OUTVCF $vcfHeader;
  while(my $x=$vcfObj->next_data_hash()){
    my $posId=$$x{CHROM}.':'.$$x{POS};
    my $pos=$$x{POS};
    my $numAlts=scalar(@{$$x{ALT}});

    # I just don't think some fields even belong here
    for(qw(ADP HET NC HOM WT)){
      delete($$x{INFO}{$_});
    }

    # Convert to haploid
    for($$x{REF}, @{$$x{ALT}}){
      $_=substr($_,0,1);
    }

    # If this is looking at a single sample, then just say that
    # the whole site passes. Then, modify the FT tag to say
    # whether it really passes.
    # The filter field will be reevaluated later in the code
    # anyway.
    #if($type eq 'single'){
      #$$x{FILTER}=['PASS'];
    #}

    # Some format tag modification:
    #   1. Make everything haploid
    #   2. Add the FT tag
    while(my($samplename,$gtypesHash)=each(%{$$x{gtypes}})){
      # Convert alleles to haploid, eg A/A => A
      $$x{gtypes}{$samplename}{GT}=substr($$x{gtypes}{$samplename}{GT},0,1);
      # Convert FILTER to FT tag
      if(!$$gtypesHash{FT} || $$gtypesHash{FT} eq '.'){
        push(@{ $$x{FORMAT} }, 'FT');
        # Transfer all FILTER to FT, if requested
        if($$settings{'filter-ft'}){
          $$x{gtypes}{$samplename}{FT}=join(";",@{$$x{FILTER}});
        }
      }
    }

    #######################
    ## MASKING
    while(my($samplename,$gtypesHash)=each(%{$$x{gtypes}})){
      my $FREQ=$$x{gtypes}{$samplename}{FREQ} || '0%';
         $FREQ='0%' if($FREQ eq '.');
         $FREQ=~s/%$//;   # remove ending percentage sign
         $FREQ=$FREQ/100; # convert to decimal
      # Mask low frequency
      if($FREQ < $altFreq && $FREQ > $inverseAltFreq){
        markFailure($x,$samplename,$settings);
      }
      # Mask for low coverage.
      # Figure out what the depth is using some common depth definitions.
      $$x{gtypes}{$samplename}{DP}||="";
      $$x{gtypes}{$samplename}{DP}||=$$x{gtypes}{$samplename}{SDP} || $$x{gtypes}{$samplename}{ADP} || 0;
      $$x{gtypes}{$samplename}{DP}=0 if($$x{gtypes}{$samplename}{DP} eq '.');
      if($$x{gtypes}{$samplename}{DP} < $$settings{min_coverage}){
        for(0 .. $numAlts-1){
          $$x{gtypes}{$samplename}{FT}='FAIL';
        }
      }
    }
    # END MASKING
    ########################

    # Reevaluate whether this whole site passes.
      # Count how many samples fail at this site.
      my $numFail=0;
      while(my($samplename,$gtypesHash)=each(%{$$x{gtypes}})){
        $numFail++ if($$gtypesHash{FT} ne 'PASS');
      }
      # If they all fail, then the site fails.
      if($numFail>=$numSamples){
        $$x{FILTER}=['FAIL'];  # Mark a failure
        $$x{ALT}=['N'];        # Mask the base
      }

    # Done! Print the modified line.
    print OUTVCF $vcfObj->format_line($x);
  }
  close OUTVCF;

  return $v;
}

sub markFailure{
  my($x,$samplename,$settings)=@_;

  # Mark which alt is the number pertaining to N
  my $nOrdinal;
  for(my $i=0;$i<@{$$x{ALT}};$i++){
    if(uc($$x{ALT}) eq 'N'){
      $nOrdinal=$i;
      last;
    }
  }
  # Add "N" if not defined since we have to mark
  # something as a failure.
  if(!defined($nOrdinal)){
    push(@{$$x{ALT}},'N');
    $nOrdinal=scalar(@{$$x{ALT}}) - 1;
  }

  # Fail a site if requested
  if($$settings{'fail-sites'}){
    $$x{FILTER}='FAIL';
    $$x{ALT}=['N'];
    # also mark each GT as the masked base
  }
  # Fail a sample at this site if requested
  if($$settings{'fail-samples'}){
    # Mark this sample as failing.
    $$x{gtypes}{$samplename}{FT}='FAIL';
    # Set the gtype equal to N.
    $$x{gtypes}{$samplename}{GT}=$nOrdinal;
  }
  if($$settings{'add-N'}){
    
  }
  if($$settings{'remove-unused-alt'}){
    
  }
}

sub usage{
  "Fixes a given VCF by adding, removing, or editing fields. This script
  will also reevaluate sites on whether they have passed filters using the FT tag.
  However, the FILTER field will simply be set to 'PASS'.
  Usage: $0 file.vcf.gz > fixed.vcf
  --numcpus      1
  --min_coverage 0
  --min_alt_frac 0
  --rename-sample ''     Rename the sample name to something. Assumes
                         that there is a single sample.
  --filter-ft            Before a site is evaluated, transfer values from
                         FILTER to the FT tag for each sample.
  --fail-sites           If a site fails, add FAIL under the FEATURE column
  --fail-samples         If a site fails, add FAIL under the FT tag for
                         the sample(s)
  --remove-unused-alt    After a site has been evaluated, if an ALT's
                         value is no sample's GT, remove it.
  "
  #--add-N                If a site or a sample fails, add N as an ALT
}

