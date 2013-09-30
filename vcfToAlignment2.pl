#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use File::Basename;
use threads;
use Thread::Queue;

my $time=time;
sub logmsg{my $duration=time-$time; $|++;print STDERR "$duration\t@_\n";$|--;}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help outfile=s reference=s badPosition=s coverage=i numcpus=i));
  $$settings{outfile}||="$0.out.fasta";
  $$settings{coverage}||=10;
  $$settings{numcpus}||=1;
  die usage($settings) if($$settings{help} || @ARGV<2);
  my $reference=$$settings{reference} or die "ERROR: need reference file\n".usage($settings);

  # organize vcf/bam files to their actual extension/type
  my(@VCF,@BAM);
  while(my $file=shift(@ARGV)){
    my($dir,$name,$ext)=fileparse($file,qw(.vcf .bam));
    if($ext=~/vcf/i){
      push(@VCF,$file);
    } elsif($ext=~/bam/i){
      push(@BAM,$file);
    } else {
      die "ERROR: Could not understand what kind of file $file is";
    }
  }

  logmsg "Finding depths for all genomes";
  my %depth=depths(\@BAM,$settings);
  logmsg "Putting all bases into a hash of arrays";
  my $refBase=findReferenceBases($reference,$settings);

  printFasta(\@VCF,\@BAM,$refBase,\%depth,$settings);
  return 0;
}

sub printFasta{
  my($VCF,$BAM,$refBase,$depthHash,$settings)=@_;

  # find all bad positions
  my (%badPos);
  if((my $file=$$settings{badPosition}) && -f $$settings{badPosition}){
    open(BAD,"<",$file) or warn "WARNING: Could not open $file for reading:$!\n No positions will be marked as 'bad'\n";
    while(<BAD>){
      chomp;
      $badPos{$_}=1;
    }
    close BAD;
  }

  # print each fasta
  my @contig=keys(%$refBase);
  for my $vcf(@$VCF){
    # find the correct bam file to check out depths
    my $genome=basename($vcf,'.vcf');
    my @bam=grep(/$genome\b/,@$BAM);
    die "ERROR: more than one bam file matches the vcf name $vcf:\n".Dumper(@bam) if(@bam>1);
    my $bam=shift(@bam);
    my $depth=$$depthHash{$bam};

    # read the vcf and make base calls
    logmsg "Printing $vcf...";
    print ">$vcf\n";
    my $vcfHash=readVcf($vcf,$settings);
    # TODO get the depth here, so that each assembly is one "unit," making multithreading easier
    for my $contig(@contig){
      my $currentPos=1;
      my @unsortedPos=keys(%{$$vcfHash{$contig}});
      for my $pos(sort{$a<=>$b} @unsortedPos){
        # fill in the sequence with reference bases until we get to the variant call
        for my $i ($currentPos .. $pos-1){
          my $posKey=$contig.'_'.$i;
          if($badPos{$posKey} || !$$depth{$posKey} || $$depth{$posKey}<$$settings{coverage}){
            print 'N'; next;
          }
          print $$refBase{$contig}[$i];
        }

        # print the alternate base, if it is above the right coverage level
        $$depth{$contig.'_'.$pos}||=0;
        if($$depth{$contig.'_'.$pos}<$$settings{coverage}){
          print 'N';
        } else {
          print $$vcfHash{$contig}{$pos};
        }

        # reset the position to the one after the alternate
        $currentPos=$pos+1;
      }

      # Use the reference genome to fill in any remaining positions after the last variant
      for my $i ($currentPos .. scalar(@{$$refBase{$contig}})-1){
        my $posKey=$contig.'_'.$i;
        if($badPos{$posKey} || !$$depth{$posKey} || $$depth{$posKey}<$$settings{coverage}){
          print 'N'; next;
        }
        print $$refBase{$contig}[$i];
      }
    }
    print "\n";
  }

  return 1;
}

sub depths{
  my($bam,$settings)=@_;
  my @thr;
  my %depth;
  my $Q=Thread::Queue->new(@$bam);
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(sub{
      my($Q,$settings)=@_;
      my $depth;
      while(defined(my $b=$Q->dequeue)){
        $$depth{$b}=covDepth($b,$settings);
      }
      return $depth;
    },$Q,$settings);
  }
  $Q->enqueue(undef) for(@thr);
  for(@thr){
    my $tmp=$_->join;
    next if(!$tmp);
    %depth=(%depth,%$tmp);
  }
  return %depth; 
}

sub covDepth{
  my($bam,$settings)=@_;
  my $depthFile="$bam.depth";
  logmsg "Finding depth for $bam";
  my %depth;
  system("samtools depth '$bam' > '$depthFile'") if(!-e $depthFile);
  die if $?;
  open(IN,$depthFile) or die "Could not open $bam in samtools:$!";
  while(<IN>){
    my($rseq,$pos,$depth)=split(/\t/);
    chomp($depth);
    $depth{$rseq."_".$pos}=$depth;
  }
  close IN;
  return \%depth;
}

sub findReferenceBases{
  my($reference,$settings)=@_;
  my $base={};
  my $in=Bio::SeqIO->new(-file=>$reference);
  while(my $seq=$in->next_seq){
    logmsg "Getting bases for contig ".$seq->id;
    # $$base{$seq->id}=[undef,split(//,$seq->seq)];
    my @seq=split(//,lc($seq->seq));
    unshift(@seq,undef);
    $$base{$seq->id}=\@seq;
  }
  return $base;
}

sub readVcf{
  my($vcf,$settings)=@_;
  my $vcfHash={};
  open(VCF,"<",$vcf) or die "ERROR: could not open vcf file $vcf:$!";
  while(<VCF>){
    next if(/^#/);
    chomp;
    my($contig,$pos,undef,undef,$alt)=split /\t/;
    $$vcfHash{$contig}{$pos}=$alt;
  }
  return $vcfHash;
}

sub usage{
  $0 = basename $0;
  "Creates an alignment of SNPs, given a set of VCFs. Positions not described in the VCF will either have N or will assume the reference nt.
  usage: $0 *.bam *.vcf -r reference.fasta -b bad.txt -c 10 > alignment.fasta
    -b bad.txt: all positions that should not be used, in format of contig_pos
    -n numcpus (default: 1)
    -r reference.fasta The reference file that all reads are already mapped against
    -c 10 the level of coverage required for a base call. An 'N' will be used if an nt is below this coverage.
  "
}
