#!/usr/bin/env perl

# TODO remove closely clustered SNPs

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
  my $settings={};
  GetOptions($settings,qw(help ambiguities min_coverage=i tempdir=s allowed=i positions=s bcfOutput=s)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{min_coverage}||=10;
  $$settings{tempdir}||="tmp";
  $$settings{allowed}||=0;
  $$settings{bcfOutput}||="$$settings{tempdir}/bcfquery";
  my $vcf=$ARGV[0];

  # read in the file with bcftools query
  my %matrix;
  my $bcfqueryFile=bcftoolsQuery($vcf,$settings);

  # Remove clustered SNPs, etc
  my @queryMatrix=filterSites($bcfqueryFile,$settings);
  
  my $header=shift(@queryMatrix);
  $header=~s/^#\s*//; # remove the hash in front of the header
  my @header=split(/\t/,$header);
  $_=~s/\[\d+\]// for(@header); # remove [number] notations for the headers
  $_=~s/:GT$//    for(@header); # remove :GT after each genotype field

  # genome names
  my @genome=@header[3..@header-1];

  for(@queryMatrix){
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
      die "ERROR: could not find the genotype field (GT)!\n". Dumper [$genome,\%F] if(!$GT);
      my $nt=$GT;
      my ($nt1,$nt2)=split(/[\/\|]/,$nt);
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
  if($$settings{positions}){
    logmsg "Sorted positions pertaining to the final (and possibly filtered) alignment are found in $$settings{positions}";
    open(POS,">",$$settings{positions}) or die "ERROR: Could not open positions file $$settings{positions}: $!";
    print POS join("\n",@POS)."\n";
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
  my $bcfquery="bcftools query -i '%TYPE=\"snp\" && %MIN(DP)>=$$settings{min_coverage}' -f '%CHROM\\t%POS\\t%REF\\t[%TGT\\t]\\n' --print-header";

  # handle the input file depending on how it is compressed or not
  my $streamInCommand="";
  if($file=~/\.gz$/){
    $streamInCommand="gunzip -c $file";
  } else {
    $streamInCommand="cat $file";
  }

  my $bcfMatrix=`$streamInCommand | $bcfquery`;
  die "ERROR running command: $!\n  $bcfquery < $file" if $?;

  # Write the bcf query to a file
  # TODO return the file instead of a large string
  open(BCFQUERY,">",$$settings{bcfOutput}) or die "ERROR: could not open $$settings{bcfOutput} for writing: $!";
  print BCFQUERY $bcfMatrix;
  close BCFQUERY;

  return $$settings{bcfOutput};
    
  #return split(/\n/,$bcfMatrix) if(wantarray);
  #return $bcfMatrix;
}

sub filterSites{
  my($bcfqueryFile,$settings)=@_;


  my ($newMatrix,$currentPos,$currentChr);
  open(BCFQUERY,"<",$bcfqueryFile) or die "ERROR: could not open bcftools query file $bcfqueryFile for reading: $!";
  while(my $bcfMatrixLine=<BCFQUERY>){
    if($bcfMatrixLine=~/^#/){
      push(@$newMatrix,$bcfMatrixLine);
      next;
    }

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
      next;
    }

    # If the SNP is far enough away, then accept it into the new matrix
    if($POS - $currentPos >= $$settings{allowed}){
      push(@$newMatrix,$bcfMatrixLine);
    }

    # update the position
    $currentPos=$POS;
  }

  return @$newMatrix if wantarray;
  return $newMatrix;
}

sub usage{
  "Multiple VCF format to alignment
  $0: reads a pooled vcf (from bcftools merge) and creates a multiple sequence alignment file.
  Usage: 
    $0 pooled.fastq.gz > aln.fasta
  --ambiguities       to allow for ambiguity letter codes other than 'N'
  --allowed 0         How close SNPs can be from each other before being thrown out
  --min_coverage 10   Minimum coverage per site before it can be considered
  --tempdir tmp       temporary directory
  --positions pos.tsv Output file with position information
  --bcfOutput out.bcf Output file with the bcftools query matrix
  "
}
