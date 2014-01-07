#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::Perl;
use File::Basename;
use threads;
use Thread::Semaphore;
use Thread::Queue;

$0=basename $0;
my $diskIoStick=Thread::Semaphore->new();
sub logmsg{print STDERR "$0: ".(caller(1))[3].": @_\n"}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help outfile=s reference=s coverage=i numcpus=i positionsFile=s));
  $$settings{outfile}||="$0.out.fasta";
  $$settings{coverage}||=10;
  $$settings{numcpus}||=1;
  die usage($settings) if($$settings{help} || @ARGV<2);
  my $reference=$$settings{reference} or die "ERROR: need reference file\n".usage($settings);

  # organize vcf/bam files to their actual extension/type
  my(@VCF,@BAM);
  while(my $file=shift(@ARGV)){
    my($dir,$name,$ext)=fileparse($file,qw(.unfiltered.vcf .vcf .bam));
    if($ext=~/vcf/i){
      push(@VCF,$file);
    } elsif($ext=~/bam/i){
      push(@BAM,$file);
    } else {
      die "ERROR: Could not understand what kind of file $file is";
    }
  }
  

  logmsg "Discovering all defined VCF positions";
  my $posArr=getVcfPositions(\@VCF,$settings);
  logmsg "Reading the reference genome";
  my $refBase=findReferenceBases($$settings{reference},$settings);
  if($$settings{positionsFile}){
    my $posThread=threads->new(\&printPositions,$posArr,$$settings{positionsFile},$settings);
    $posThread->detach; # this will be done soon after running, so it should be done before the script is done
  }

  # kick off VCF=>fasta threads
  my $Q=Thread::Queue->new;
  my $printQ=Thread::Queue->new;
  my @thr;
  $thr[$_]=threads->new(\&vcfToFastaWorker,\@BAM,$posArr,$refBase,$Q,$printQ,$settings) for(0..$$settings{numcpus}-1);
  $Q->enqueue(@VCF);

  # create printer worker thread
  my $printer=threads->new(sub{
    my($Q,$settings)=@_;
    while(defined(my $thing=$Q->dequeue)){
      print $thing;
    }
  },$printQ,$settings);

  # wrap up threads
  $Q->enqueue(undef) for(@thr);
  $_->join for(@thr);
  # wrap up printer thread
  $printQ->enqueue(undef);
  $printer->join;

  return 0;
}

sub vcfToFastaWorker{
  my($BAM,$posArr,$refBase,$Q,$printQ,$settings)=@_;

  while(defined(my $vcf=$Q->dequeue)){
    my $genome=basename($vcf,qw(.unfiltered.vcf .vcf)); # done for each vcf...
    # get the correct bam for the vcf
    my @bam=grep(/$genome\b/,@$BAM);
    die "ERROR: there are many bams that fit the description $genome: ".join(" ",@bam) if(@bam>1);
    my $bam=$bam[0];

    # get the depth of this bam
    my $coverage=covDepth($bam,$settings);

    logmsg "Printing the fasta entry for $vcf";
    my $fasta=fastaEntry($vcf,$posArr,$coverage,$refBase,$settings);
    $printQ->enqueue($fasta);
  }

  return 0;
}

# Reads all VCFs to find a union of all positions.
# Returns a sorted array of positions, with each position in the form of contig_pos
sub getVcfPositions{
  my($VCF,$settings)=@_;
  my %pos;
  my %bad;
  for my $file(@$VCF){
    my ($vcf,$bad)=readVcf($file,$settings);
    while(my($contig,$posHash)=each(%$vcf)){
      for my $pos(keys(%$posHash)){
        $pos{$contig}{$pos}=1;
      }
    }
    # Have to keep track of and not delete at this step, because all VCFs need to be read first.
    while(my($contig,$posHash)=each(%$bad)){
      for my $pos(keys(%$posHash)){
        $bad{$contig}{$pos}=1;
      }
    }
  }

  # Find where there are any 'bad' sites and remove them
  while(my($contig,$posHash)=each(%bad)){
    for my $pos(keys(%$posHash)){
      delete($pos{$contig}{$pos}) if($pos{$contig}{$pos});
    }
  }

  # put the hash to a sorted array
  my (@pos);
  while(my($contig,$posIndex)=each(%pos)){
    for my $pos(keys(%$posIndex)){
      push(@pos,$contig.'_'.$pos);
    }
  }
  my @sortedPos=sort {
    my ($contigA,$posA)=split(/_/,$a);
    my ($contigB,$posB)=split(/_/,$b);
    return $contigA cmp $contigB if($contigA ne $contigB);
    return $posA <=> $posB;
  } @pos;
    
  return \@sortedPos;
}

sub readVcf{
  my($vcf,$settings)=@_;
  my $vcfHash={};
  my $badHash={};
  $diskIoStick->down;
  open(VCF,"<",$vcf) or die "ERROR: could not open vcf file $vcf:$!";
  while(<VCF>){
    next if(/^#/);
    chomp;
    my($contig,$pos,undef,$ref,$alt)=split /\t/;
    $$vcfHash{$contig}{$pos}=$alt;

    if($ref eq '*' || $alt eq '*' || length($ref)>1 || length($alt)>1){
      $$badHash{$contig}{$pos}=1;
    }
  }
  close VCF;
  $diskIoStick->up;
  return ($vcfHash,$badHash) if wantarray;
  return $vcfHash;
}

sub fastaEntry{
  my ($vcf,$posArr,$coverage,$refBase,$settings)=@_;
  my $v=readVcf($vcf,$settings); # The tradeoff for using lower memory: reading the vcf file a second time here
  my $fasta="";
  # print defline
  $fasta.= ">$vcf\n";
  for my $posKey(@$posArr){
    my ($contig,$pos)=split(/_/,$posKey);
    my $nt=$$v{$contig}{$pos} || 'N';

    # If the base caller didn't say anything about this base, see if there is 
    # enough coverage to call it the reference base.
    if($nt eq 'N'){
      if($$coverage{$posKey} && $$coverage{$posKey} >= $$settings{coverage}){
        $nt=$$refBase{$contig}[$pos];
      }
    }
    $fasta.=$nt;
  }
  $fasta.="\n";
  return $fasta;
}

sub covDepth{
  my($bam,$settings)=@_;
  my $depthFile="$bam.depth";
  my %depth;
  #system("samtools depth '$bam' > '$depthFile'") if(!-e $depthFile);
  die if $?;
  open(IN,"samtools depth '$bam' | ") or die "Could not open $bam in samtools:$!";
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

sub printPositions{
  my($pos,$file,$settings)=@_;
  # make a copy of positions so that it doesn't mess up the rest of the script
  my @pos=@$pos;
  open(POS,">",$file) or die "ERROR: Could not write to positions file $file: $!";
  for(@pos){
    print POS $_."\n";
  }
  close POS;
  logmsg "Sorted positions have been printed to $file";
  return 1;
}

sub usage{
  "Creates an alignment of SNPs, given a set of VCFs
  usage: $0 *.bam *.vcf -o alignment.fasta -r reference.fasta
    -n numcpus (default: 1)
    -coverage 10 The minimum coverage allowed to accept the reference base, if the base caller didn't call a position.
    -p positions.txt To output positional information to this file. Each line of the positions file corresponds to the respective position in the fasta alignment file.
  "
}
