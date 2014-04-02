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
  GetOptions($settings,qw(help outfile=s reference=s coverage=i numcpus=i positionsFile=s table=s minimumDistance=i)) or die;
  $$settings{outfile}||="$0.out.fasta";
  $$settings{coverage}||=10;
  $$settings{numcpus}||=1;
  $$settings{minimumDistance}||=0;
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
  my $tableQ=Thread::Queue->new; # for holding snp table results (if requested)
  my @thr;
  $thr[$_]=threads->new(\&vcfToFastaWorker,\@BAM,$posArr,$refBase,$Q,$printQ,$tableQ,$settings) for(0..$$settings{numcpus}-1);
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

  printSnpTable($posArr,$tableQ,$$settings{table},$settings) if($$settings{table});

  return 0;
}

sub printSnpTable{
  my($posArr,$tableQ,$file,$settings)=@_;
  $tableQ->enqueue(undef);
  open(TSV,">$file") or die "ERROR: Could not open table file $file:$!";
  print TSV join("\t",@$posArr)."\n";    
  while(defined(my $row=$tableQ->dequeue)){
    print TSV $row;
  }
  close TSV;
}

sub vcfToFastaWorker{
  my($BAM,$posArr,$refBase,$Q,$printQ,$tableQ,$settings)=@_;

  while(defined(my $vcf=$Q->dequeue)){
    my $genome=basename($vcf,qw(.unfiltered.vcf .vcf)); # done for each vcf...
    # get the correct bam for the vcf
    my @bam=grep(/$genome\b/,@$BAM);
    die "ERROR: there are many bams that fit the description $genome: ".join(" ",@bam) if(@bam>1);
    my $bam=$bam[0];

    # get the depth of this bam
    my $coverage=covDepth($bam,$settings);

    logmsg "Printing the fasta entry for $vcf";
    my ($fasta,$tableRow)=fastaEntry($vcf,$posArr,$coverage,$refBase,$settings);
    $printQ->enqueue($fasta);
    $tableQ->enqueue($tableRow);
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
      push(@pos,$contig.'_'.$pos);
    }
  }
  my @sortedPos=sort {
    my ($contigA,$posA)=split(/_/,$a);
    my ($contigB,$posB)=split(/_/,$b);
    return $contigA cmp $contigB if($contigA ne $contigB);
    return $posA <=> $posB;
  } @pos;

  # return the positions if not removing by min distance
  return \@sortedPos if(!$$settings{minimumDistance});

  # if there is a min distance specified, remove those SNPs that are too close
  my @newSortedPos;
  my $numPos=@sortedPos;
  my ($currentDistance,$currentContig,$currentPos)=(0,split(/_/,$sortedPos[0]));
  for(my $i=1;$i<$numPos;$i++){
    my($contig,$pos)=split(/_/,$sortedPos[$i]);
    my($prevContig,$prevPos)=split(/_/,$sortedPos[$i-1]);
    if($contig ne $currentContig){
      $currentContig=$contig;
      next;
    }
    my $distance=$pos - $prevPos;
    # keep distantly-related positions
    if($distance >= $$settings{minimumDistance}){
      push(@newSortedPos,$sortedPos[$i]);
    }
  }

  return \@newSortedPos;
}

sub readVcf{
  my($vcf,$settings)=@_;
  my $vcfHash={};
  $diskIoStick->down; # mark that one process is using the disk
  open(VCF,"<",$vcf) or die "ERROR: could not open vcf file $vcf:$!";
  while(<VCF>){
    next if(/^#/);
    chomp;
    my($contig,$pos,undef,$ref,$alt)=split /\t/;
    $$vcfHash{$contig}{$pos}=$alt;

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

sub fastaEntry{
  my ($vcf,$posArr,$coverage,$refBase,$settings)=@_;
  my $v=readVcf($vcf,$settings); # The tradeoff for using lower memory: reading the vcf file a second time here
  my $fasta="";
  my $table="";
  # print defline
  $fasta.= ">$vcf\n";
  $table.="$vcf\t" if($$settings{table});
  for my $posKey(@$posArr){
    my ($contig,$pos)=split(/_/,$posKey);
    my $nt=$$v{$contig}{$pos} || 'N';

    # If the base caller didn't say anything about this base, see if there is 
    # enough coverage to call it the reference base.
    # Do not check lowercase N because those are not due to low coverage.
    if($nt eq 'N'){
      if($$coverage{$posKey} && $$coverage{$posKey} >= $$settings{coverage}){
        $nt=$$refBase{$contig}[$pos];
      }
    }
    $fasta.=$nt;

    if($$settings{table}){
      $table.="$nt\t";
    }
  }
  $fasta.="\n";
  $table.="\n" if($$settings{table});
  return ($fasta,$table) if wantarray;
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
    die "ERROR: the reference sequence contig has an underscore in it. You must map and call snps from contigs without underscores or else there will be internal problems in this script. Problem was with contig ".$seq->id if($seq->id =~/_/);
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
  "Creates an alignment of SNPs, given a set of VCFs. Output is in fasta format.
  usage: $0 *.bam *.vcf -r reference.fasta > alignment.fasta
    Note: multiple nucleotide polymorphic sites and indels are changed to 'n'. Use filterVcf.pl to help retain some of these positions.
    -n numcpus (default: 1)
    -coverage 10 The minimum coverage allowed to accept the snp. Or, the reference base, if the base caller didn't call a position.
    -p positions.txt To output positional information to this file. Each line of the positions file corresponds to the respective position in the fasta alignment file.
    -t table.txt To output SNP calls to a table
    -m 0 The minimum distance allowed between two hqSNPs. Those that are too close will be all-together removed from the MSA.
  "
}
