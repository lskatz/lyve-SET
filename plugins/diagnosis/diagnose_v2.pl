#!/usr/bin/env perl
# Authors: Khushbu Patel and Lori Gladney 3-16-18
# Lyve-SET diagnosis/QC script
# Edited, Lee Katz 2018-05-08
# Usage: ./hqSNP_diagnose.pl projectdir /path/to/reference/ref.fasta

use warnings;
use strict;
use Bio::SeqIO;
use Data::Dumper qw/Dumper/;
use Getopt::Long qw/GetOptions/;
use File::Basename qw/basename dirname/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../../lib";
use LyveSET qw/logmsg/;

# Bring in all script dependencies here
$ENV{PATH} = "$ENV{PATH}:$RealBin/../../scripts";

# Global variables
my @maskingHeaders=qw(alignment ambiguous_base_count ambiguous_base_percent GC% contig_length);
my @bamHeaders=qw(bam total_reads mapped_reads mapped_reads_ge30 unmapped_reads percent_mapped percent_mapped_gt30 proper_pairs inward_pairs outward_pairs pairs_with_other_orientation percentage_proper_pairs insert_size_avg number_sites_with_coverage_lt20 breadth_coverage total_sites_covered_ge20 avg_coverage_per_site perc_reference_coverage);
my @vcfHeaders=qw(vcf PASS AF0.95 AF0.95;str10 DP20 AF0.95;DP20 DP20;AF0.95 DP20;AF0.95;str10 AF0.95;DP20;str10 AF0.95;isIndel AF0.95;RF0.95 AF0.95;isIndel;RF0.95 isIndel isIndel;AF0.95 RF0.95 isIndel;RF0.95 str10 str10;AF0.95 str10;DP20;AF0.95 indelError);

# Run the main subroutine, and the exit code for the script
# is the return integer from main().
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  die usage() if($$settings{help});
  my $project = $ARGV[0];
  die "ERROR: need project directory\n".usage() if(!$project);
  $project=~s/\/+$//; # trim off any leading slashes for readability

  #.....Setting the paths to fetch files!.........#
  my $path_cleaned_reads = "$project/reads";
  my $vcf_path           = "$project/vcf";
  my $bam_path           = "$project/bam";
  my $out_aln_path       = "$project/msa/out.aln.fasta";
  my $out_info_path      = "$project/msa/out.informative.fasta";
  my $reference          = "$ARGV[1]";
  my $out_pool_vcf		 = "$project/msa/out.pooled.vcf.gz";
  

  if(!-e $reference){
    die "ERROR: could not find reference file $reference";
  }

  # the reference size
  logmsg "Reading the reference genome";
  my $total_ref_bases = gsize($reference);
  
  # out.pool.vcf info
  print "# Out.pool.vcf Info\n";
  #system("zcat $out_pool_vcf | grep -v '^#' | cut -f 7| sort | uniq -c");
  my $poolVCF = glob("$project/msa/out.pooled.vcf.gz");
  print join("\t",@vcfHeaders)."\n";
  
  # set defaults
  my %filter_code;
  for(@vcfHeaders){
    $filter_code{$_}//=".";
	}
	
  open(pool_vcf,"zcat $out_pool_vcf | ") or die "ERROR: could not zcat $out_pool_vcf: $!";
  while(<pool_vcf>){
    next if(/^#/);
    chomp;

    my(@F) = split('\t');
	$filter_code{$F[6]}++;
	}
	
	close pool_vcf;
	

  print $filter_code{$_}."\t" for(@vcfHeaders);
  print "\n\n";
	
  # Masking information on both pseudo alignments
  logmsg "Reading pseudo alignments";
  my $maskingInfoHash = masking($out_info_path);
  my $maskingAlnHash  = masking($out_aln_path);

  # Print the masking report
  print "# Full alignment\n";
  print join("\t",@maskingHeaders)."\n";
  for my $m (@$maskingAlnHash) {
    print join("\t", $$m{id}, $$m{countN}, $$m{percentN}, $$m{percentGC}, $$m{length})."\n";
  }
  print "\n";
  print "# Informative alignment\n";
  print join("\t",@maskingHeaders)."\n";
  for my $m (@$maskingInfoHash) {
    print join("\t", $$m{id}, $$m{countN}, $$m{percentN}, $$m{percentGC}, $$m{length})."\n";
  }

  # Print the bam report
  print "# BAM report\n";
  print join("\t",@bamHeaders)."\n";
  foreach my $bam(glob("$bam_path/*.bam")){
    logmsg "Processing bam $bam...";
    my $bamInfo = bamInfo($bam,$path_cleaned_reads,$total_ref_bases,$vcf_path);
    print $$bamInfo{$_}."\t" for(@bamHeaders);
    print "\n";
  }

  # VCF report
  print "# VCF report\n";
  print join("\t",@vcfHeaders)."\n";
  foreach my $vcf(glob("$vcf_path/*.vcf.gz")){
    my $vcfInfo = vcfInfo($vcf,$path_cleaned_reads,$total_ref_bases,$vcf_path);
    print $$vcfInfo{$_}."\t" for(@vcfHeaders);
    print "\n";
  }

  return 0;
}

sub bamInfo{
  my($bam,$path_cleaned_reads,$total_ref_bases,$vcf_path)=@_;

  # Default return hash
  my %bamInfo=(bam=>$bam);
  for my $field (@bamHeaders){
    $bamInfo{$field}//=".";
  }

  # Example filename: bam/sample1.fastq.gz-reference.sorted.bam
  # Example basename: sample1
  my $basename = basename($bam,qw(.fastq.gz-reference.sorted.bam));	
  $bamInfo{bam}=$basename;
  

  #Fetching total number of mapped reads
  #NOTE: saving system() output to perl variable, returns 0
  #it sets the exit status to perl variable; Hence use ``
  my $total_mapped_reads = `samtools view -c -F 4 $bam` + 0; # cast to int; remove newline
  die if $?;
  #Fetching number of mapped reads with quality greater than or equal to 30
  my $mapped_reads_ge_30 = `samtools view -q 30 $bam | wc -l` + 0;
  die if $?;
  #Number of reads < Q30
  my $mapped_reads_lt_30 = $total_mapped_reads - $mapped_reads_ge_30;

  $bamInfo{mapped_reads}     =$total_mapped_reads + 0;
  $bamInfo{mapped_reads_ge30}=$mapped_reads_ge_30 + 0;
  $bamInfo{mapped_reads_lt30}=$mapped_reads_lt_30 + 0;

  #Calculating total number of reads from clean_reads files
  # Gather hopefully just one file pertaining to this bam
  my $temp = $basename =~ s/.fastq*.+//;
  my @file = glob("$path_cleaned_reads"."/"."$basename*gz");
  if(@file > 1){
    logmsg "Warning: found more than one fastq file for $bam";
  } elsif(@file < 1){
    logmsg "Warning: could not find a fastq file for $bam";
    return \%bamInfo;
  }
  my $file = join(" ",@file);
  #Gettting total number of lines in fastq files: the
  #number of lines divided by 4.
  my $line_count= `zcat $file | wc -l`;
  die if $?;
  my $tot_reads = $line_count/4;
  $bamInfo{total_reads}=$tot_reads;

  $bamInfo{unmapped_reads} = $tot_reads - $total_mapped_reads; # Calculating number of reads that are not mapped
  $bamInfo{percent_mapped} = ($total_mapped_reads/$tot_reads) * 100; #Calculating percentage of mapped reads
  $bamInfo{percent_mapped_gt30} = ($mapped_reads_ge_30/$total_mapped_reads) * 100; #Calculating percentage of mapped reads > Q30

  #my @bamHeaders=qw(bam total_reads mapped_reads mapped_reads_gt30 unmapped_reads percent_mapped percent_mapped_gt30 proper_pairs inward_pairs outward_pairs pairs_with_other_orientation percentage_proper_pairs insert_size_avg number_sites_with_coverage_lt20 breadth_coverage total_sites_covered_ge20 avg_coverage_per_site);

  # Run samtools one time and fill in some more variables
  my $stats_line_counter=0;
  open(SAMTOOLS_STATS, "samtools stats $bam | ") or die "ERROR running samtools stats $bam: $!";
  while(<SAMTOOLS_STATS>){
    last if($stats_line_counter++ > 40);
    chomp;

    if     (/reads properly paired:\s*(\d+)/){
      $bamInfo{proper_pairs}=$1;
      $bamInfo{percentage_proper_pairs}=sprintf("%0.2f",$bamInfo{proper_pairs}/$bamInfo{mapped_reads_ge30} * 100);
    } elsif(/inward oriented pairs:\s*(\d+)/){
      $bamInfo{inward_pairs}=$1;
    } elsif(/outward oriented pairs:\s*(\d+)/){
      $bamInfo{outward_pairs}=$1;
    } elsif(/pairs with other orientation:\s*(\d+)/){
      $bamInfo{pairs_with_other_orientation}=$1;
    } elsif(/insert size average:\s*(\d+)/){
      $bamInfo{insert_size_avg}=$1;
    }
  }
  close SAMTOOLS_STATS;

  # Percentage of coverage of the reference genome with depth >=20
  my $sites_lt20=0;
  my $sites_ge20=0;
  my $coverage_count=0;
  my $bam_sites=0;
  open(SAMTOOLS_DEPTH, "samtools depth $bam |") or die "ERROR running samtools depth $bam: $!";
  while(<SAMTOOLS_DEPTH>){
    chomp;
    my($chr,$pos,$depth)=split(/\t/,$_);
    $coverage_count+=$depth;
    $bam_sites++;
    if($depth >= 0 && $depth < 20){
      $sites_lt20++;
    } else {
      $sites_ge20++;
    }
  }
  close SAMTOOLS_DEPTH;
  $bamInfo{number_sites_with_coverage_lt20}=$sites_lt20;
  $bamInfo{total_sites_covered_ge20}=$sites_ge20;
  $bamInfo{breadth_coverage}=$bam_sites;
  $bamInfo{avg_coverage_per_site}=sprintf("%0.2f",$coverage_count/$total_ref_bases);
  $bamInfo{perc_reference_coverage}= ((($total_ref_bases - $sites_lt20)/$total_ref_bases)*100);


  return \%bamInfo;
}

sub vcfInfo{
  my($vcf,$path_cleaned_reads,$total_ref_bases,$vcf_path)=@_;
  my $basename = basename($vcf,qw(-reference.vcf.gz));


=cut
  print("Filters applied for sites that do not pass -
  a. DP20 = Depth is less than 20, the user-set coverage threshold
  b. RF0.95 = Reference variant consensus is less than 0.95, the user-set threshold
  c. AF 0.95 = Allele variant consensus is less than 0.95, the user-set threshold
  d. isIndel = Indels are not used for analysis in Lyve-SET
  e. masked = This site was masked using a bed file or other means
  f. str10 = Less than 10% or more than 90% of variant supporting reads on one strand
  g. indelError = Likey artifcat due to indel reads at this position\n\n");
=cut

  # set defaults
  my %filter_code = (vcf=>$basename);
  for(@vcfHeaders){
    $filter_code{$_}//=".";
	
  }
  

  open(VCF,"zcat $vcf | ") or die "ERROR: could not zcat $vcf: $!";
  while(<VCF>){
    next if(/^#/);
    chomp;

    my(@F) = split('\t');
	$filter_code{$F[6]}++;
	
  }
  close VCF;

  return \%filter_code;
} # end of vcfInfo subroutine


#subroutine to calculate masking from out.aln and out.informative file
sub masking() {
  my ($file) = @_;

  #create one SeqIO object to read in
  my $seq_in = Bio::SeqIO->new(
               -file   => "<$file",
               -format => "fasta",
               );


  my @unsorted_id;  #array to hold sequence IDs and associated info

  #Save each entry in the input file.
  #seq_in is the object that holds the file
  #next seq passes each contig to the sequence object $seq
  while (my $seq = $seq_in->next_seq()) { # Loops while condition is true

    #Initialize variables; start the counter at 0 for each contig that comes through the loop

    my $id= $seq->id; #Get the header Id of each contig
    my $length=$seq->length;  #Declare variable for finding the length of each contig sequence

    my $sequence   = $seq->seq;   #Each contig seq is treated as a string
    my $countN     = $sequence =~ s/(N)//sgi; # number of ambiguous bases (N)
       $countN   ||=0;  # set to zero if anything
    my $percentN   = sprintf("%0.2f", ($countN / $length) *100);
    my $countGC    = $sequence =~ s/[GC]//sgi;
    my $percentGC  = sprintf("%0.2f",$countGC/$length * 100);

    my $full_id= "$id %Ambiguous_bases(N)= $percentN Ambiguous_base_count(N)= $countN %GC= $percentGC Contig_length= $length";

    #Push each scalar into the unsorted array (identifiers, sequence and percent ambiguous bases)
    push (@unsorted_id,{
        id         => $id,
        full_id    => $full_id,
        length     => $length,
        countN     => $countN,
        percentN   => $percentN,
        percentGC  => $percentGC,
      }
    );
  }

  #Loop through the array.. Print the contig sequences in order of GC content (lowest to highest)
  my @sorted_id = sort{ $$a{percentGC} <=> $$b{percentGC} } @unsorted_id;
  return \@sorted_id;

} #end of masking subroutine

#Get total number of bases in the Reference
sub gsize(){
  my($reference)=@_;
  #create one SeqIO object to read in
  my $seq_in = Bio::SeqIO->new(
               -file   => "<$reference",
               -format => 'fasta',
               );



  my $id;        #scalar to hold sequence identifier
  my $string;    #scalar to hold the sequence
  my $length =0;

  #write each entry in the input file to STDOUT, while the condition is true
  #seq_in is the object that holds the file
  #next seq passes each contig to the sequence object $seq
  while (my $seq = $seq_in->next_seq ()) {
    #Initialize variables; start the counter at 0 for each contig that comes through the loop
    $length +=$seq->length;  #Declare variable for finding the length of each contig sequence
  }

  return $length;

}

sub usage{
  "Runs a diagnosis on a finished Lyve-SET run
  Usage: ./hqSNP_diagnose.pl projectdir reference.fasta
  "
}
