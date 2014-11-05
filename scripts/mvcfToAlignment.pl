#!/usr/bin/env perl

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
  GetOptions($settings,qw(help ambiguities));
  die usage() if($$settings{help});
  
  my $header=<>;
  chomp($header);
  $header=~s/^#\s*//; # remove the hash in front of the header
  my @header=split(/\t/,$header);
  $_=~s/\[\d+\]// for(@header); # remove [number] notations for the headers
  $_=~s/:GT$//    for(@header); # remove :GT after each genotype field

  # genome names
  my @genome=@header[3..@header-1];

  # read the matrix into memory before printing the aln
  my %matrix;
  while(<>){
    chomp;
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

sub usage{
  "Multiple VCF format to alignment
  $0: reads a bcftools query output and creates a multiple sequence alignment file. Required columns: GT (genotype), REF, CHROM, POS
  Usage: bcftools query [...] | $0 > aln.fasta
  --ambiguities  to allow for ambiguity letter codes other than 'N'
  "
}
