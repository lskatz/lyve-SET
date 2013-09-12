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
sub logmsg{my $duration=time-$time; $|++;print "$duration\t@_\n";$|--;}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help outfile=s reference=s badPosition=s coverage=i allowedFlanking=i numcpus=i));
  $$settings{outfile}||="$0.out.fasta";
  $$settings{coverage}||=10;
  $$settings{numcpus}||=1;
  $$settings{allowedFlanking}||=0;
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

  logmsg "Ok! Done getting depths and reference sequence information. Converting vcf to fasta alignment now";
  my ($fastaStr,$pos)=vcfToFasta(\@VCF,\@BAM,$refBase,\%depth,$settings);

  # print fasta to a file
  open(FASTA,">",$$settings{outfile}) or die "Cannot open file for writing: $$settings{outfile}: $!";
  print FASTA $fastaStr; 
  close FASTA;
  print "Alignment is in $$settings{outfile}\n";

  # print which positions for each nt in the alignment to a file
  open(POS,">","$$settings{outfile}.pos.txt") or die "Cannot open file for writing: $$settings{outfile}.pos.txt:$!";
  print POS join("\t",@$pos)."\n";
  close POS;
  print "List of variant positions is in $$settings{outfile}.pos.txt\n";

  return 0;
}

sub vcfToFasta{
  my($VCF,$BAM,$refBase,$depth,$settings)=@_;

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

  # find all positions that the VCFs encompass
  # filter out all bad positions, defined by badPosition
  my (%pos,@genome);
  for my $vcf(@$VCF){
    open(VCF,"<",$vcf) or die "Could not open $vcf for reading:$!";
    my $genome=basename($vcf,qw(.vcf));
    push(@genome,$genome);
    while(<VCF>){
      next if(/^\s*$/ || /^#/);
      chomp;
      my @F=split /\t/;
      my $posKey=join("_",@F[0..1]);
      die if(!$posKey);
      next if($badPos{$posKey});
      my $alt=$F[4];
      $pos{$posKey}{$genome}=$alt;
    }
    close VCF;
  }
  
  ## sort the positions on two fields
  my @pos=sort({
    my ($cB,$posB)=split(/_/,$b);
    my ($cA,$posA)=split(/_/,$a);
    $cA cmp $cB || 
    $posA <=>$posB;
  } keys(%pos));

  # for every possible position, see what each genome's variant call is
  for my $genome(@genome){
    my @bam=grep(/$genome\b/,@$BAM);
    die "ERROR: there are many bams that fit the description $genome: ".join(" ",@bam) if(@bam>1);
    my $bam=$bam[0];
    $$depth{$genome}=$$depth{$bam} if(!$$depth{$genome});
    my %genomeDepth=%{$$depth{$genome}};
    for my $posKey(@pos){
      next if($pos{$posKey}{$genome} || $badPos{$posKey});
      my($rseq,$pos)=split(/_/,$posKey);
      if($genomeDepth{$posKey} && $genomeDepth{$posKey} >= $$settings{coverage}){
        $pos{$posKey}{$genome}=$$refBase{$rseq}[$pos];
      } else {
        $pos{$posKey}{$genome}="N";

        # Better yet, don't delete and allow for the complete picture
        #$badPos{$posKey}=1;
        #delete($pos{$posKey}); 
      }
    }
  }

  # now that we have base calls for each genome at each site, make a fasta
  my $allowedFlanking=$$settings{allowedFlanking};
  my ($fasta,@legitPos);
  my $numPositions=@pos;
  for my $genome(@genome){
    my ($prevContig,$prevPos);
    $fasta.=">$genome\n";
    for(my $i=0;$i<$numPositions;$i++){
      my $posKey=$pos[$i];
      next if($badPos{$posKey});

      # do not accept positions that are too close together
      my($contig,$position)=split /_/,$posKey;
      next if($prevPos && ($contig eq $prevContig) && ($position-$allowedFlanking < $prevPos));

      $fasta.=$pos{$posKey}{$genome};
      push(@legitPos,$posKey);
      $prevContig=$contig;
      $prevPos=$position;
    }
    $fasta.="\n";
  }
  return ($fasta,\@legitPos);
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
  print "Finding depth for $bam\n";
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
    logmsg $seq->id;
    my @seq=split(//,$seq->seq);
    unshift(@seq,undef);
    $$base{$seq->id}=\@seq;
  }
  return $base;
}

sub usage{
  "Creates an alignment of SNPs, given a set of VCFs
  usage: $0 *.bam *.vcf -o alignment.fasta -r reference.fasta -b bad.txt
    -b bad.txt: all positions that should not be used, in format of contig_pos
    -a allowed flanking in bp (default: 0)
      nucleotides downstream of another snp this many bp away will not be accepted
    -n numcpus (default: 1)
  "
}
