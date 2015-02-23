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
  my $settings={};
  GetOptions($settings,qw(help tempdir=s)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{tempdir}||="tmp";

  my($matrix)=@ARGV;

  bcfqueryToFasta($matrix,$settings);

  return 0;
}

sub bcfqueryToFasta{
  my($bcfqueryFile,$settings)=@_;
  
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
  close BCFQUERY;
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
  $out->close;

  return 1;
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

sub usage{
  "Multiple VCF format to alignment
  $0: reads a tab-delimted file created by 'bcftools query'.
  Usage: 
    $0 bcfmatrix.tsv > alignment.fasta
    $0 < bcfmatrix.tsv > alignment.fasta # read stdin
  --tempdir tmp       temporary directory
  "
}
