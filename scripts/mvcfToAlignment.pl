#!/usr/bin/env perl

# TODO put the bcftools query command into this script so that it actually reads from an multiple-genome vcf

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
  GetOptions($settings,qw(help ambiguities min_coverage=i tempdir=s));
  die usage() if($$settings{help} || !@ARGV);
  $$settings{min_coverage}||=10;
  $$settings{tempdir}||="tmp";
  my $vcf=$ARGV[0];

  # read in the file with bcftools query
  my %matrix;
  my @queryMatrix=bcftoolsQuery($vcf,$settings);
  chomp(@queryMatrix);
  
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
  my $bcfMatrix=`$bcfquery < $file `;
  return split(/\n/,$bcfMatrix) if(wantarray);
  return $bcfMatrix;
}

sub usage{
  "Multiple VCF format to alignment
  $0: reads a pooled vcf (from bcftools merge) and creates a multiple sequence alignment file.
  Usage: 
    bcftools *.vcf.gz -O z > pooled.fastq.gz # prerequisite
    $0 pooled.fastq.gz > aln.fasta           # $0 command
  --ambiguities  to allow for ambiguity letter codes other than 'N'
  --min_coverage 10 Minimum coverage per site before it can be considered
  "
}
