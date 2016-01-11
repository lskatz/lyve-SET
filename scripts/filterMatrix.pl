#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use File::Basename;
use File::Temp qw/tempdir/;
use List::Util qw/min max/;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/logmsg/;
use lib "$FindBin::RealBin/../lib/lib/perl5";
use Number::Range;

$0=fileparse $0;

exit main();
sub main{
  my $settings={ambiguities=>1,invariant=>1, 'invariant-loose'=>1};
  GetOptions($settings,qw(help ambiguities! invariant! tempdir=s allowed|allowedFlanking=i mask=s@ invariant-loose!)) or die $!;
  die usage() if($$settings{help});
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{allowed}||=0;
  $$settings{mask}||=[];
  $$settings{invariant}=0 if(!$$settings{'invariant-loose'});

  my($in)=@ARGV;
  $in||=""; # avoid undefined warnings

  # Explain the settings for this script.
  my $settingsString="";
  for (sort {$a cmp $b} keys(%$settings)){
    my $str=$$settings{$_};
    if(ref($str) eq 'ARRAY'){
      $str=join("\t",@$str);
    }
    $settingsString.=join("\t=\t",$_,$str)."\n";
  }
  logmsg "Filtering and other settings for this script:\n$settingsString";
  filterSites($in,$settings);

  return 0;
}

# Filter the BCF query file into a new one, if any filters were given
sub filterSites{
  my($bcfqueryFile,$settings)=@_;

  # Open the unfiltered BCF query file
  my $fp;
  if($bcfqueryFile){
    open($fp,"<",$bcfqueryFile) or die "ERROR: could not open bcftools query file $bcfqueryFile for reading: $!";
  } else {
    $fp=*STDIN;
  }

  # If there are any coordinates explicitly listed for masking,
  # combine them all into a single range object, one per seqname.
  my $maskedRanges=readBedFiles($bcfqueryFile,$$settings{mask},$settings);

  # Read in the header with genome labels
  my $header=<$fp>;
  print $header;  # Print the header right away so that it can be saved correctly before it's changed
  $header=~s/^\s+|^#|\s+$//g; # trim and remove pound sign
  my @header=split /\t/, $header;

  # Read the actual content with variant calls
  my ($currentPos,$currentChr);
  while(my $bcfMatrixLine=<$fp>){
    my $lineLength=length($bcfMatrixLine);  # for seeking
    chomp $bcfMatrixLine;

    # start by assuming that this is a high-quality site
    my $hqSite=1;

    # get the fields from the matrix
    my($CONTIG,$POS,$REF,@GT)=split(/\t/,$bcfMatrixLine);
    my $numAlts=@GT;

    for(my $i=0;$i<$numAlts;$i++){
      $GT[$i]=diploidGtToHaploid($GT[$i],$REF,$settings);
    }

    # get the values of the contig/pos if they aren't set
    if(!defined($currentChr) || $CONTIG ne $currentChr){
      # If the contig's value has switched, then undefine the position and start comparing on the new contig
      undef($currentChr);
      undef($currentPos);

      # start anew with these coordinates
      $currentChr=$CONTIG;
      $currentPos=$POS;
      seek($fp,-$lineLength,1);
      next;
    }
    # Alter the bcf line after the seeking line, so that the correct seek length can be established.
    $bcfMatrixLine=join("\t",$CONTIG,$POS,$REF,@GT);

    # Mask any site found in the BED files
    $hqSite=0 if(defined($$maskedRanges{$CONTIG}) && $$maskedRanges{$CONTIG}->inrange($POS));

    # High-quality sites are far enough away from each other, as defined by the user.
    # This step also assumes that there is a SNP right before the contig starts
    # and does not allow a SNP to start at the beginning of the contig if --allowed
    # is turned on.
    # TODO: have this pseudo SNP at the end of the contig too.
    $hqSite=0 if($POS - $currentPos < $$settings{allowed});

    # Simply get rid of any site that consists of all Ns
    # TODO make this optional somehow.
    my $is_allNs=1;
    # TODO test this code with grep
    #   $is_allNs=0 if(grep(!/N/,@GT));
    for(my $i=0;$i<$numAlts;$i++){
      if($GT[$i]!~/N/i){
        $is_allNs=0;
        last; # save a nanosecond
      }
    }
    $hqSite=0 if($is_allNs);

    # The user can specify that high quality sites are those where every site is defined
    # (ie through --noambiguities)
    if(!$$settings{ambiguities}){
      for(my $i=0;$i<$numAlts;$i++){
        if($GT[$i]!~/[ATCG]/i){
          $hqSite=0;
          last;
        }
      }
    }

    # Remove any invariant site, if specified.
    # invariant means "keep invariant sites"
    # !invariant means "remove invariant sites"
    # However, there is no point in looking for these
    # sites if all alts are masked as Ns.
    if(!$is_allNs && !$$settings{invariant}){

      # Need a ref base that will be in the MSA that
      # is not ambiguous for decent comparison.
      # The bases in the MSA are only coming from alts.
      my $altRefIndex=0;
      my $altRef='N';
      while($altRef=~/N/i){
        $altRef=$GT[$altRefIndex];
        die "Internal error: could not find unambiguous base at this site: $CONTIG:$POS" if($altRefIndex > 99999);
        $altRefIndex++;
      }
      die "Internal error: alt reference not found at $CONTIG:$POS" if(!$altRef);

      # Start off assuming that the site is not variant,
      # until proven otherwise.
      my $is_variant=0;
      for(my $i=1;$i<$numAlts;$i++){
        if($GT[$i] ne $altRef){
          next if(!$$settings{'invariant-loose'} &&  $GT[$i]=~/N/i);
          $is_variant=1;
          last; # save a nanosecond of time here
        }
      }
      $hqSite=0 if(!$is_variant);
    }
    
    # print out the results
    print $bcfMatrixLine ."\n" if($hqSite);

    # update the position
    $currentPos=$POS;
  }
  close $fp;

  return 1;
}

sub diploidGtToHaploid{
  my($gt,$REF,$settings)=@_;

  my($gt1,$gt2)=split(/\//,$gt);
  $gt2||=$gt1;
  for($gt1,$gt2){
    $_=$REF if($_ eq ".");
  }
  if($gt1 ne $gt2){
    $gt="N";
  } else {
    $gt=$gt1;
  }

  return $gt;
}

# Turn bed-defined ranges into range objects
sub readBedFiles{
  my($bcfqueryFile,$maskFile,$settings)=@_;
  my %range;

  # Find what seqnames there are out there first,
  # so that every expected range object will be defined.
  if(-e $bcfqueryFile){
    my $command="grep -v '^#' $bcfqueryFile | cut -f 1 | sort | uniq";
    my $seqname=`$command`;
    die "ERROR getting sequence names from $bcfqueryFile: $!\n  $command" if $?;
    my @seqname=split(/\n/,$seqname);
    chomp(@seqname);
    undef($seqname);

    # Initialize ranges
    $range{$_}=Number::Range->new() for(@seqname);
  }

  # Read the bed files
  for my $bed(@$maskFile){
    open(BED,$bed) or die "ERROR: could not open $bed: $!";
    while(<BED>){
      chomp;
      my($seqname,$start,$stop)=split(/\t/);
      $range{$seqname}=Number::Range->new() if(!$range{$seqname});
      my $lo=min($start,$stop);
      my $hi=max($start,$stop);
      $range{$seqname}->addrange("$lo..$hi");
    }
    close BED;
  }

  return \%range;
}

sub usage{
  "Multiple VCF format to alignment
  $0: filters a bcftools query matrix. The first three columns of the matrix are contig/pos/ref, and the next columns are all GT.
  Usage: 
    $0 bcftools.tsv > filtered.tsv
    $0 < bcftools.tsv > filtered.tsv # read stdin
  --noambiguities               Remove any site with an ambiguity
                                (i.e., complete deletion)
  --noinvariant                 Remove any site whose alts are all equal.
  --noinvariant-loose           Invokes --noinvariant, but does not consider 
                                ambiguities when deciding whether a site is 
                                variant. Sites can consist of only REF and N.
  --allowed           0         How close SNPs can be from each other before 
                                being thrown out
  --tempdir           tmp       temporary directory
  --mask              file.bed  BED-formatted file of regions to exclude.
                                Multiple --mask flags are allowed for multiple
                                bed files.
  "
}
