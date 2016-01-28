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
  my @thresholdOpts=qw(min_coverage|min-coverage|coverage=i min_alt_frac|min-alt-frac=s DF|min-DF=i DR|min-DR=i DP4=i);
  my @lyvesetOpts=qw(help tempdir=s numcpus=i);
  GetOptions($settings,qw(fail-all fail-samples fail-sites rename-sample=s type=s pass-until-fail),@thresholdOpts,@lyvesetOpts) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1,CLEANUP=>1);
  mkdir $$settings{tempdir};

  # Set zero as the min threshold if not specified.
  $$settings{$_}||=0 for(qw(min_coverage min_alt_frac DF DR DP4));
  # Some parameter checking
  die "ERROR: min_alt_frac must be between 0 and 1, inclusive" if($$settings{min_alt_frac} < 0 || $$settings{min_alt_frac} > 1);
  $$settings{type} && logmsg "WARNING: --type is no longer used in Lyve-SET >= v1.3";
  
  # Shorthand options
  if($$settings{'fail-all'}){
    $$settings{'fail-samples'}=1;
    $$settings{'fail-sites'}=1;
  }
  if(my $DP4=$$settings{'DP4'}){
    $$settings{$_}=$DP4 for(qw(DF DR));
  }

  my $VCF=$ARGV[0];
  die usage($settings) if(!$VCF || $$settings{help});
  
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
  my $reevaluated=reevaluateSites($uncompressed,$settings);
  logmsg "Done fixing it up!";

  return $reevaluated;
}

sub addStandardTags{
  my($compressedVcf,$settings)=@_;
  my $v="$$settings{tempdir}/filledInAnAcRef.vcf";
  logmsg "Excluding indels from the analysis first...";
  
  system("bcftools view $compressedVcf --exclude-types indels > $v");
  die "ERROR: could not run bcftools on $compressedVcf." if $?;

  # rename the sample if requested
  if(my $newname=$$settings{'rename-sample'}){
    logmsg "Renaming the sample to $newname";
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
  $vcfObj->add_header_line({key=>'FORMAT', ID=>'DP4', Number=>'4', Type=>'Integer', Description=>"Number of forward ref alleles; reverse ref; forward non-ref; reverse non-ref alleles."}) if(!@{ $vcfObj->get_header_line(key=>'FORMAT',ID=>'DP4') });
  for my $tag(qw(RDF RDR ADF ADR)){
    $vcfObj->add_header_line({key=>'FORMAT', ID=>$tag, Number=>1, Type=>'Integer', Description=>"See: DP4"}) if(!@{ $vcfObj->get_header_line(key=>'FORMAT',ID=>$tag) });
  }

  # start printing
  my $vcfHeader=$vcfObj->format_header();
  print OUTVCF $vcfHeader;
  my $posCounter=0; # How many positions, not the genomic position
  while(my $x=$vcfObj->next_data_hash()){
    # Add in Lyve-SET goodies
    $x=fixVcfLine($x,$settings);

    $posCounter++;
    my $posId=$$x{CHROM}.':'.$$x{POS};
    my $pos=$$x{POS};
    my $numAlts=scalar(@{$$x{ALT}});
    my @samplename=keys(%{$$x{gtypes}});
    my $numSamples=@samplename;

    #######################
    ## MASKING
    for my $samplename(@samplename){
      # Figure out an alt-frequency between 0 and 1
      my $FREQ=$$x{gtypes}{$samplename}{FREQ} || '0%';
         $FREQ='0%' if($FREQ eq '.');
         $FREQ=~s/%$//;   # remove ending percentage sign
         $FREQ=$FREQ/100; # convert to decimal
      # Mask low frequency
      if($FREQ < $altFreq && $FREQ > $inverseAltFreq){
        if($$settings{'fail-samples'}){
          $$x{gtypes}{$samplename}{GT}=$$x{altIndex}{N};
          $$x{gtypes}{$samplename}{FT}='FAIL';
        }
        if($$settings{'fail-sites'}){
          $$x{gtypes}{$samplename}{GT}=$$x{altIndex}{N};
          $$x{FILTER}=['FAIL'];
        }
      }

      # Mask for low coverage.
      # Figure out what the depth is using some common depth definitions.
      $$x{gtypes}{$samplename}{DP}||="";
      $$x{gtypes}{$samplename}{DP}||=$$x{gtypes}{$samplename}{SDP} || $$x{gtypes}{$samplename}{ADP} || 0;
      $$x{gtypes}{$samplename}{DP}=0 if($$x{gtypes}{$samplename}{DP} eq '.');
      if($$x{gtypes}{$samplename}{DP} < $$settings{min_coverage}){
        if($$settings{'fail-samples'}){
          $$x{gtypes}{$samplename}{GT}=$$x{altIndex}{N};
          $$x{gtypes}{$samplename}{FT}='FAIL';
        }
        if($$settings{'fail-sites'}){
          $$x{gtypes}{$samplename}{GT}=$$x{altIndex}{N};
          $$x{FILTER}=['FAIL'];
        }
      }

      # Mask for DP4
      my($RDF, $RDR, $ADF, $ADR)=split(/\s*,\s*/,$$x{gtypes}{$samplename}{DP4});
      for($RDF, $RDR, $ADF, $ADR){
        $_=0 if($_ eq '.');
      }
      if(($RDF + $ADF) < $$settings{DF} || ($RDR + $ADR) < $$settings{DR}){
        if($$settings{'fail-samples'}){
          $$x{gtypes}{$samplename}{GT}=$$x{altIndex}{N};
          $$x{gtypes}{$samplename}{FT}='FAIL';
        }
        if($$settings{'fail-sites'}){
          $$x{gtypes}{$samplename}{GT}=$$x{altIndex}{N};
          $$x{FILTER}=['FAIL'];
        }
      }

    }
    # END MASKING
    ########################

    # Reevaluate whether this whole site passes.
    # Count how many samples fail at this site.
    my $numFail=0;
    for my $samplename(@samplename){
      $numFail++ if($$x{gtypes}{$samplename}{FT} ne 'PASS');
    }
    # If they all fail, then the site fails.
    if($numFail==$numSamples){
      $$x{FILTER}=['FAIL'];  # Mark a failure
      # Setting ALT equal to ['N'] is unnessessary because N is already
      # present and all other ALTs will be eliminated when they are 
      # seen as unused.
    }

    # Remove alternate bases that are not used.
    # This includes a fake REF ALT that this script
    # introduces and a fake N that may or may not
    # be used depending on site/sample filters.
    removeUnusedAlts($x,$settings);

    # Done! Print the modified line.
    print OUTVCF $vcfObj->format_line($x);
    
    if($posCounter % reportEvery == 0){
      logmsg "Looked at $posCounter positions so far";
    }
  }
  close OUTVCF;

  return $v;
}

# Remove alternate bases that are not used.
# If there is more than one alt, and a dot
# is one of the ALT values, remove it.
sub removeUnusedAlts{
  my($x,$settings)=@_;

  my @samplename=keys(%{$$x{gtypes}});

  # Mark which ALTs are in common with GTs.
  # Keys are integer GT values; values are 1 to mark true.
  my %usedAlt;
  # Explicitly mark the reference base as being used.
  # It will be taken into account later.
  $usedAlt{0}=1;
  for my $samplename(@samplename){
    $usedAlt{$$x{gtypes}{$samplename}{GT}}=1;
  }
  
  # For each unused ALT, delete it from the VCF position.
  # This has to be done in reverse-ordinal position such that
  # the splice doesn't trip over itself.
  my $numAlts=@{$$x{ALT}};
  for(my $ordinal=$numAlts-1;$ordinal>=0;$ordinal--){
    next if($usedAlt{$ordinal});

    # Since we know that at this point in the loop, the ALT is
    # not used, delete it and fix the GT values greater
    # than or equal to this ordinal
    #logmsg "Splicing $ordinal ($$x{POS}: $$x{ALT}[$ordinal]) on ".Dumper $$x{ALT};
    splice(@{$$x{ALT}},$ordinal,1);
    # Reassign GT integers: it should be subtracted if it is at 
    # least as large as the ALT ordinal that was just deleted.
    for my $samplename(@samplename){
      $$x{gtypes}{$samplename}{GT}-- if($$x{gtypes}{$samplename}{GT} > $ordinal);
    }
  }

  # Remove the REF ALT (index 0)
  shift(@{$$x{ALT}});
  # Add back the REF ALT as a dot if there are no other ALTs
  unshift(@{$$x{ALT}},'.') if(@{$$x{ALT}} == 0);

 
  #die Dumper $x;
  #die Dumper $x if($$x{POS}==550);
}

# Add in the reference nt and the ambiguous "N" nt explicitly.
# Add additional information such as ALT indexing.
sub fixVcfLine{
  my($x,$settings)=@_;

  # Make sure I don't mess with $x by reference.
  my %y=%$x;
  # Initial declarations
  my @samplename=keys(%{$y{gtypes}});
  my $numSamples=@samplename;
  
  # Make the ordinals correct: add the reference base as
  # an "alt" so the index matches up in @ALT.
  # The REF's GT value is 0, so this is an unshift function.
  if(!defined($y{ALT}[0]) ||
    ($y{ALT}[0] ne '.' && $y{ALT}[0] ne $y{REF})
  ){
    unshift(@{$y{ALT}},$y{REF});
  }

  # Change $y by reference if there are any 
  # '.' in the ALT or GT fields
  for(my $i=0;$i<@{$y{ALT}};$i++){
    # Make the ref base explicit
    $y{ALT}[$i]=$y{REF} if($y{ALT}[$i] eq '.');
    # Convert to haploid and uppercase
    $y{ALT}[$i]=uc(substr($y{ALT}[$i],0,1));
  }
  $y{REF}=uc($y{REF});

  # Innocent until proven guilty
  if($$settings{'pass-until-fail'}){
    $y{FILTER}=['PASS'];
  } else {
    # If we are not assuming innocent until guilty,
    # and if there is no value set yet, then it is FAIL
    if(!defined($y{FILTER}[0]) || $y{FILTER}[0] eq '.'){
      $y{FILTER}=['FAIL'];
    }
  }

  # Index the FORMAT tags
  $y{formatHash}={};
  $y{formatHash}{$_}=1 for(@{$y{FORMAT}});

  # Fix up the samples' tags
  for my $samplename(@samplename){
    # Convert alleles to haploid first, eg 1/1 => 1
    $y{gtypes}{$samplename}{GT}=substr($y{gtypes}{$samplename}{GT},0,1);
    # I have noticed a dot instead of an integer in the GT field.
    # This needs to be replaced with an ordinal representing
    # the REF which is 0.
    $y{gtypes}{$samplename}{GT}="0" if($y{gtypes}{$samplename}{GT} eq '.');

    # Add an FT tag to describe whether the sample's GT passes.
    if(!$y{gtypes}{$samplename}{FT}){
      # Only add it to this array if it hasn't been added yet.
      # We don't want to add it once per sample.
      if(!$y{formatHash}{FT}++){
        push(@{ $y{FORMAT} }, 'FT');
      }
      # must set a default value in case it doesn't get set
      $y{gtypes}{$samplename}{FT}='FAIL';
    }
    # add RDF RDR ADF ADR tags
    for my $tag(qw(RDF RDR ADF ADR)){
      next if(defined($y{gtypes}{$samplename}{$tag}));
      if(!$y{formatHash}{$tag}++){
        push(@{ $y{FORMAT} },$tag);
      }
      $y{gtypes}{$samplename}{$tag}='.';
    }

    # Add DP4 tag
    if(!defined($y{gtypes}{$samplename}{DP4})){
      if(!$y{formatHash}{DP4}++){
        push(@{ $y{FORMAT} }, 'DP4');
      }
      # Take the DP4 value from INFO if possible
      if($y{INFO}{DP4}){
        $y{gtypes}{$samplename}{DP4}=$y{INFO}{DP4};
        delete($y{INFO}{DP4});
      } 
      # If it wasn't defined before, set it up
      $y{gtypes}{$samplename}{DP4}=join(",",
        $y{gtypes}{$samplename}{RDF},
        $y{gtypes}{$samplename}{RDR},
        $y{gtypes}{$samplename}{ADF},
        $y{gtypes}{$samplename}{ADR}
      );

    }
    # Convert FILTER to FT tag
    if($y{gtypes}{$samplename}{FT} eq '.'){
      # Transfer all FILTER to FT, if there is no value set yet
      $y{gtypes}{$samplename}{FT}=join(";",@{$y{FILTER}});
    }
    # Innocent until proven guilty
    if($$settings{'pass-until-fail'}){
      $y{gtypes}{$samplename}{FT}='PASS';
    }else{
      # If we are not assuming innocent until guilty,
      # and if there is no value set yet, then it is FAIL
      if($y{gtypes}{$samplename}{FT} eq '.'){
        $y{gtypes}{$samplename}{FT}='FAIL';
      }
    }
  }
  
  # Mark which alt is the number pertaining to N
  # or any other nucleotide.
  # REF is the zeroth ALT in VCF format and in this array
  # as it is set up in this subroutine.
  my $ordinal={};
  my $numAlts=@{$y{ALT}};
  for(my $i=0;$i<$numAlts;$i++){
    $$ordinal{$y{ALT}[$i]}=$i;
  }

  # Add "N" if not defined since we have to mark
  # something as a failure for this subroutine.
  if(!$$ordinal{N}){
    $$ordinal{N}=$numAlts;
    push(@{$y{ALT}},'N');
  }
  $y{altIndex}=$ordinal;

  # Fix up the info column.
  # I just don't think some fields even belong here
  for(qw(ADP HET NC HOM WT DP4)){
    delete($y{INFO}{$_});
  }
    
  #die Dumper \%y if($y{POS}==1);
  #die Dumper \%y if($y{POS}==550);

  return \%y;
}

sub usage{
  my($settings)=@_;
  $settings||={};
  my $usage="Fixes a given VCF by adding, removing, or editing fields. 
  Without any options, useless INFO tags and ALT nts are
  removed.  Additionally, each site will be failed if all samples
  have failed in the FT tag field
  Usage: $0 file.vcf.gz > fixed.vcf
  --min_coverage 0
  --min_alt_frac 0
  --DP4          0       Shorthand for --min-DF x --min-DR x where x 
                         is an integer. If these individual options are 
                         not found in the VCF but DP4 is, the actual 
                         DP4 value will be used.

  --fail-all             Short for '--fail-sites --fail-samples'
  --fail-sites           If a site fails, add FAIL under the FEATURE column
  --fail-samples         If a site fails, add FAIL under the FT tag for
                         the sample(s).
  --pass-until-fail      The 'innocent-until-proven guilty' argument.
                         Assume every site and sample passes unless
                         proven otherwise. Without this argument, previous
                         FAIL filters will remain in place.

  --help                 for more help
  ";

  return $usage if(!$$settings{help});

  $usage.="
  --numcpus      1       Num CPUS currently has no effect on this script.
  --tempdir      /tmp/   Choose a temporary directory (optional)
  --rename-sample ''     Rename the sample name to something. Assumes
                         that there is a single sample.

  --min-DF       0       Minimum number of supporting bases on FWD strand
  --min-DR       0       Minimum number of supporting bases on REV strand
  ";
  return $usage;
}

