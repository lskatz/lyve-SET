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
my $diskIoStick;
sub logmsg{print STDERR "$0: ".(caller(1))[3].": @_\n"}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help reference=s coverage=i numcpus=i)) or die;
  $$settings{coverage}||=10;
  $$settings{numcpus}||=1;
  die usage($settings) if($$settings{help} || @ARGV<2);
  my $reference=$$settings{reference} or die "ERROR: need reference file\n".usage($settings);

  # organize vcf/bam files to their actual extension/type
  my(@VCF,@BAM);
  while(my $file=shift(@ARGV)){
    my($dir,$name,$ext)=fileparse($file,qw(.unfiltered.vcf .vcf .vcf.gz .bam));
    if($ext=~/vcf$/i || $ext=~/vcf\.gz$/){
      push(@VCF,$file);
    } elsif($ext=~/bam$/i){
      push(@BAM,$file);
    } else {
      die "ERROR: Could not understand what kind of file $file is";
    }
  }

  # Mark how many files we can read from the drive at once
  $diskIoStick=Thread::Semaphore->new($$settings{numcpus});

  logmsg "Reading the reference genome";
  my $refBase=findReferenceBases($$settings{reference},$settings);
  logmsg "Discovering all defined VCF positions";
  my $posArr=getVcfPositions(\@VCF,$settings);

  # kick off VCF=>matrix threads
  my $Q=Thread::Queue->new;
  # Start off the matrix with a positions header
  my $printQ=Thread::Queue->new(vcfHeader($VCF[0],$settings));
  my @thr;
  $thr[$_]=threads->new(\&vcfToTableWorker,\@BAM,$posArr,$refBase,$Q,$printQ,$settings) for(0..$$settings{numcpus}-1);
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

sub vcfToTableWorker{
  my($BAM,$posArr,$refBase,$Q,$printQ,$settings)=@_;

  while(defined(my $vcf=$Q->dequeue)){
    my $genome=basename($vcf,qw(.unfiltered.vcf .vcf .unfiltered.vcf.gz .vcf.gz)); # done for each vcf...
    # get the correct bam for the vcf
    my @bam=grep(/$genome\b/,@$BAM);
    die "ERROR: there are many bams that fit the description $genome: ".join(" ",@bam) if(@bam>1);
    my $bam=$bam[0];
    if(!$bam || !-e $bam){
      logmsg "ERROR: I could not find a bam file to match against your vcf file $vcf but was expecting something with the prefix $genome";
      logmsg "These are the possible bam files that you have given me to choose from: ".Dumper $BAM;
      die;
    }

    # get the depth of this bam
    my $coverage=covDepth($bam,$settings);

    logmsg "Printing the entries for $vcf";
    my ($tableRows)=snpsEntry($vcf,$posArr,$coverage,$refBase,$settings);
    $printQ->enqueue($tableRows);
  }

  return 0;
}

# Reads all VCFs to find a union of all positions.
# Returns a sorted array of positions, with each position in the form of contig_pos
sub getVcfPositions{
  my($VCF,$settings)=@_;
  my %pos;
  for my $file(@$VCF){
    my $vcf=readVcf($file,$settings);
    while(my($contig,$posHash)=each(%$vcf)){
      for my $pos(keys(%$posHash)){
        $pos{$contig}{$pos}=1;
      }
    }
  }

  # put the hash to a sorted array
  my (@pos);
  while(my($contig,$posIndex)=each(%pos)){
    for my $pos(keys(%$posIndex)){
      push(@pos,$contig.':'.$pos);
    }
  }
  my @sortedPos=sort {
    my ($contigA,$posA)=split(/:/,$a);
    my ($contigB,$posB)=split(/:/,$b);
    return $contigA cmp $contigB if($contigA ne $contigB);
    return $posA <=> $posB;
  } @pos;

  # return the positions if not removing by min distance
  return \@sortedPos;
}

sub vcfHeader{
  my($vcf,$settings)=@_;
  if($vcf=~/\.gz$/){ # gzipped
    open(VCF,"gunzip -c '$vcf' | ") or die "ERROR: could not open vcf file $vcf:$!";
  } else {
    open(VCF,"<",$vcf) or die "ERROR: could not open vcf file $vcf:$!";
  }
  my $header="";
  while(my $line=<VCF>){
    last if($line=~/^[^#]/); # don't go past the header lines
    # The line with only one hash is the column names.
    # Change it to the five-column format
    if($line=~/^#[^#]/){
      $line="#".join("\t",qw(VCF CHROM POS REF ALT))."\n";
    }

    $header.=$line;
  }
  close VCF;
  return $header;
}

sub readVcf{
  my($vcf,$settings)=@_;
  my $vcfHash={};
  $diskIoStick->down; # mark that one process is using the disk

  if($vcf=~/\.gz$/){ # gzipped
    open(VCF,"gunzip -c '$vcf' | ") or die "ERROR: could not open vcf file $vcf:$!";
  } else {
    open(VCF,"<",$vcf) or die "ERROR: could not open vcf file $vcf:$!";
  }
  while(<VCF>){
    next if(/^#/);
    chomp;
    my($contig,$pos,undef,$ref,$alt)=split /\t/;
    $$vcfHash{$contig}{$pos}=$alt;

    if($alt eq '.'){
      $$vcfHash{$contig}{$pos}=$ref;
    } else {
      $$vcfHash{$contig}{$pos}=$alt;
    }

    # Indels will just be an N because it is too difficult to deal with those.
    if($ref eq '*' || $alt eq '*' || length($ref)>1 || length($alt)>1){
      # lowercase N to mask it later, when trying to recover other bases for 
      # low coverage. This N is not due to low coverage.
      $$vcfHash{$contig}{$pos}='n'; 
    }
  }
  close VCF;
  $diskIoStick->up; # mark that one process is no longer using the disk
  return $vcfHash;
}

# Find all the SNPs for a given genome at hqPositions
sub snpsEntry{
  my ($vcf,$posArr,$coverage,$refBase,$settings)=@_;
  my $v=readVcf($vcf,$settings); # The tradeoff for using lower memory: reading the vcf file a second time here
  my $table="";
  # print defline
  for my $posKey(@$posArr){
    my ($contig,$pos)=split(/:/,$posKey);
    my $nt=$$v{$contig}{$pos} || 'N';

    # If the base caller didn't say anything about this base, see if there is 
    # enough coverage to call it the reference base.
    # Do not check lowercase N because those are not due to low coverage.
    if($nt eq 'N'){
      if($$coverage{$posKey} && $$coverage{$posKey} >= $$settings{coverage}){
        $nt=$$refBase{$contig}[$pos];
      }
    }

    $table.=join("\t",$vcf,$contig,$pos,$$refBase{$contig}[$pos],$nt)."\n";
  }
  return $table;
}

sub covDepth{
  my($bam,$settings)=@_;
  my $depthFile="$bam.depth";
  my %depth;
  #system("samtools depth '$bam' > '$depthFile'") if(!-e $depthFile);
  #die if $?;
  open(IN,"samtools depth '$bam' | ") or die "Could not open $bam in samtools:$!";
  while(<IN>){
    my($rseq,$pos,$depth)=split(/\t/);
    chomp($depth);
    $depth{$rseq.":".$pos}=$depth;
  }
  close IN;
  return \%depth;
}

sub findReferenceBases{
  my($reference,$settings)=@_;
  my $base={};
  my $in=Bio::SeqIO->new(-file=>$reference);
  while(my $seq=$in->next_seq){
    die "ERROR: the reference sequence contig has an underscore in it. You must map and call snps from contigs without underscores or else there will be internal problems in this script. Problem was with contig ".$seq->id if($seq->id =~/_/);
    logmsg $seq->id;
    my @seq=split(//,$seq->seq);
    unshift(@seq,undef);
    $$base{$seq->id}=\@seq;
  }
  return $base;
}


sub usage{
  "Creates a matrix of SNPs, given a set of VCFs. Output is in tall/skinny format with five columns: genome, seqname, position, ref, alt.  Output is similar to VCF format.
  usage: $0 *.bam *.vcf [*.vcf.gz] -r reference.fasta > hqSNPs.tsv
    Note: multiple nucleotide polymorphic sites and indels are changed to 'n'. Use filterVcf.pl to help retain some of these positions.
    -r file.fasta  The reference genome assembly
    -n 1           The number of cpus to use
    -c 10          The minimum coverage allowed to accept the snp. Or, the reference base, if the base caller didn't call a position.
    -h             help menu
  "
}
