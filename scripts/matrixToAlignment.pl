#!/usr/bin/env perl

# Convert a bcfquery matrix to an alignment

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use File::Basename;

$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();
sub main{
  my $settings={filter=>1};
  GetOptions($settings,qw(help ambiguities! tempdir=s allowed=i prefix=s filter!)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{tempdir}||="tmp";
  $$settings{allowed}||=0;

  # output files
  #$$settings{prefix}||="$0.out";
  #$$settings{bcfOutput}||="$$settings{prefix}.bcfquery";
  #$$settings{filteredBcf}||="$$settings{prefix}.filtered.bcfquery";
  #$$settings{fasta}   ||="$$settings{prefix}.aln.fas";
  #$$settings{positions}||="$$settings{prefix}.pos";

  #my $vcf=$ARGV[0];

  # read in the file with bcftools query
  #logmsg "Running bcftools query to find all SNPs";
  #my %matrix;
  #my $bcfqueryFile=bcftoolsQuery($vcf,$settings);

  # Remove clustered SNPs, etc, unless the user doesn't want it
  my $filteredMatrix=$bcfqueryFile;
  if($$settings{filter}){
    logmsg "Filtering for clustered SNPs or any other specified filters";
    $filteredMatrix=filterSites($bcfqueryFile,$settings);
  }

  my $fasta=bcfqueryToFasta($filteredMatrix,$settings);

  return 0;
}

sub bcfqueryToFasta{
  my($bcfqueryFile,$settings)=@_;
  my $outfile=$$settings{fasta};
  
  open(BCFQUERY,$bcfqueryFile) or die "ERROR: could not open $bcfqueryFile: $!";

  # Process the header line
  my $header=<BCFQUERY>;
  $header=~s/^#\s*//; # remove the hash in front of the header
  $header=~s/^\s+|\s+$//g; # trim whitespace
  my @header=split(/\t/,$header);
  $_=~s/\[\d+\]// for(@header); # remove [number] notations for the headers
  $_=~s/:GT$//    for(@header); # remove :GT after each genotype field

  # genome names from the header
  my($headerChr,$headerPos,$headerRef,@genome)=@header;

  # Figure out the basecall for each pos on each genome
  logmsg "Finding basecalls based on the genotype in the pooled VCF file";
  # TODO: multithread
  my (%matrix,$i);
  while(<BCFQUERY>){
    s/^\s+|\s+$//g; # trim whitespace
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
      die "ERROR: could not find the genotype field (GT) for genome '$genome'\n". Dumper [$genome,\%F] if(!$GT);
      my $nt;
      my ($nt1,$nt2)=split(/[\/\|]/,$GT); # Split on the forward slash in case it shows as diploid
      $nt2=$nt1 if(!$nt2);                # In case it is a haploid-style pooled VCF
      if($nt1 eq $nt2){
        $nt=$nt1;
      } elsif($$settings{'ambiguities'}){
        $nt=rev_iub([$nt1,$nt2],$settings);
        logmsg "Warning: GT $GT was interpreted as an ambiguity $nt at $genome/$extra{CHROM}/$extra{POS}";
      } else{
        $nt="N";
      }

      # ALT is defined as REF if it is a dot
      $nt=$extra{REF} if($nt eq '.' || $nt eq ',');
      # take care of indels or ambiguities
      $nt="N" if(length($nt)!=1 || $nt!~/^\w$/i);

      $matrix{$genome}{$extra{CHROM}}{$extra{POS}}=$nt;
    }
    $i++;
    if($i % 10000 == 0){
      logmsg "Finished calling $i SNPs";
    }
  }
  # DONE getting it into memory!

  # START output
  # get all the sorted positions
  my @CHROM=keys(%{ $matrix{$genome[0]} });
  my @POS;
  # Sort number 1: sort by chromosome
  for my $chr(sort{$a cmp $b} @CHROM){
    # Sort number 2: sort by numerical position
    push(@POS,"$chr:$_") for(sort{$a<=>$b} keys(%{ $matrix{$genome[0]}{$chr} }));
  }

  # create the alignment
  logmsg "Found all positions; converting to fasta";
  my $out=Bio::SeqIO->new(-file=>">$outfile",-format=>"fasta");
  for my $genome(@genome){
    my $seq="";
    for(@POS){
      my($chr,$pos)=split(/:/,$_);
      $seq.=$matrix{$genome}{$chr}{$pos};
    }
    my $seqObj=Bio::Seq->new(-seq=>$seq,-id=>$genome);
    $out->write_seq($seqObj);
  }
  $out->close;

  if($$settings{positions}){
    logmsg "Sorted positions pertaining to the final (and possibly filtered) alignment are found in $$settings{positions}";
    open(POS,">",$$settings{positions}) or die "ERROR: Could not open positions file $$settings{positions}: $!";
    print POS "$_\n" for(@POS);
    #print POS join("\n",@POS)."\n";
    close POS;
  }
  
  return 0;
}

sub rev_iub{
  my($ntArr,$settings)=@_;
  die "ERROR: ntArr is not an array here:\n".Dumper $ntArr if(ref($ntArr ne "ARRAY"));
  return $$ntArr[0] if(@$ntArr==1);

  # Figure out the ambiguity now
  # from bioperl, http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html
  my %REV_IUB = (
          A    => 'A',
          T    => 'T',
          U    => 'U',
          C    => 'C',
          G    => 'G',
          AC   => 'M',
          AG   => 'R',
          AT   => 'W',
          CG   => 'S',
          CT   => 'Y',
          GT   => 'K',
          ACG  => 'V',
          ACT  => 'H',
          AGT  => 'D',
          CGT  => 'B',
          ACGT => 'N',
          N    => 'N'
      );

  my $sortedNt=join("",sort{$a cmp $b} @$ntArr);
  my $ambiguity=$REV_IUB{$sortedNt};
  die "ERROR: could not understand nt for ambiguities!\n".Dumper($ntArr) if(!$ambiguity);
  return $ambiguity;
}

sub bcftoolsQuery{
  my($file,$settings)=@_;
  my $matrix="$$settings{tempdir}/$$.matrix";
  my $bcfquery="bcftools query -i '%TYPE=\"snp\"' -f '%CHROM\\t%POS\\t%REF\\t[%TGT\\t]\\n' --print-header";
  #my $bcfquery="bcftools query -i '%TYPE=\"snp\" && %MIN(DP)>=$$settings{min_coverage}' -f '%CHROM\\t%POS\\t%REF\\t[%TGT\\t]\\n' --print-header";

  # handle the input file depending on how it is compressed or not
  my $streamInCommand="";
  if($file=~/\.gz$/){
    $streamInCommand="gunzip -c $file";
  } else {
    $streamInCommand="cat $file";
  }

  # Open a bcftools query stream
  open(QUERYSTREAM,"$streamInCommand | $bcfquery |") or die "ERROR running command: $!\n  $bcfquery < $file";
  # Write the bcf query to a file
  open(BCFQUERY,">",$$settings{bcfOutput}) or die "ERROR: could not open $$settings{bcfOutput} for writing: $!";
  while(<QUERYSTREAM>){
    print BCFQUERY $_;
  }
  close BCFQUERY;
  close QUERYSTREAM;

  # return the file instead of a large string
  return $$settings{bcfOutput};
}

# Filter the BCF query file into a new one, if any filters were given
sub filterSites{
  my($bcfqueryFile,$settings)=@_;

  my $outfile=$$settings{filteredBcf};
  open(FILTEREDOUT,">",$outfile) or die "ERROR: could not open filtered output bcf query file: $!";

  open(BCFQUERY,"<",$bcfqueryFile) or die "ERROR: could not open bcftools query file $bcfqueryFile for reading: $!";

  # Read in the header with genome labels
  my $header=<BCFQUERY>;
  print FILTEREDOUT $header;  # Print the header right away so that it can be saved correctly before it's changed
  $header=~s/^\s+|^#|\s+$//g; # trim and remove pound sign
  my @header=split /\t/, $header;

  # Read the actual content with variant calls
  my ($currentPos,$currentChr);
  while(my $bcfMatrixLine=<BCFQUERY>){

    # get the fields from the matrix
    my($CONTIG,$POS,$REF,@GT)=split(/\t/,$bcfMatrixLine);
    
    # get the values of the contig/pos if they aren't set
    if(!defined($currentChr) || $CONTIG ne $currentChr){
      # If the contig's value has switched, then undefine the position and start comparing on the new contig
      undef($currentChr);
      undef($currentPos);

      # start anew with these coordinates
      $currentChr=$CONTIG;
      $currentPos=$POS;
      seek(BCFQUERY,-length($bcfMatrixLine),1);
      next;
    }

    # If the SNP is far enough away, then accept it into the new matrix
    if($POS - $currentPos >= $$settings{allowed}){
      print FILTEREDOUT $bcfMatrixLine;
    }

    # update the position
    $currentPos=$POS;
  }
  close BCFQUERY;

  close FILTEREDOUT;
  return $outfile;
}

sub usage{
  "Multiple VCF format to alignment
  $0: reads a pooled vcf (from bcftools merge) and creates a multiple sequence alignment file.
  Usage: 
    $0 pooled.fastq.gz --prefix out
  --ambiguities       to allow for ambiguity letter codes other than 'N'
  --nofilter          Do not filter out any SNPs (overrides --ambiguities)
  --allowed 0         How close SNPs can be from each other before being thrown out
  --tempdir tmp       temporary directory
  --prefix ./prefix   The prefix for all output files
  "
}
