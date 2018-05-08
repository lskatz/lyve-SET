#!/usr/bin/env perl  
# Authors: Khushbu Patel and Lori Gladney 3-16-18
# Lyve-SET diagnosis/QC script
# Edited, Lee Katz 2018-05-08

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

# Run the main subroutine, and the exit code for the script
# is the return integer from main().
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  die usage() if($$settings{help});
  my $project = $ARGV[0];
  die "ERROR: need project directory\n".usage() if(!$project);

  #.....Setting the paths to fetch files!.........#
  my $path_cleaned_reads = "$project/reads";
  my $vcf_path           = "$project/vcf";
  my $bam_path           = "$project/bam";
  my $out_aln_path       = "$project/msa/out.aln.fasta";
  my $out_info_path      = "$project/msa/out.informative.fasta";
  my $reference          = "$project/reference/reference.fasta";

  # the reference size         
  logmsg "Reading the reference genome";
  my $total_ref_bases = gsize($reference);

  # Masking information on both pseudo alignments
  logmsg "Reading pseudo alignments";
  my $maskingInfoHash = masking($out_info_path);
  my $maskingAlnHash  = masking($out_aln_path);

  # Print the masking report
  print join("\t",qw(full_alignment ambiguous_base_count ambiguous_base_percent GC% contig_length))."\n";
  for my $m (@$maskingAlnHash) {
    print join("\t", $$m{id}, $$m{countN}, $$m{percentN}, $$m{percentGC}, $$m{length})."\n";
  }
  print "\n";
  print join("\t",qw(informative_alignment ambiguous_base_count ambiguous_base_percent GC% contig_length))."\n";
  for my $m (@$maskingInfoHash) {
    print join("\t", $$m{id}, $$m{countN}, $$m{percentN}, $$m{percentGC}, $$m{length})."\n";
  }

  foreach my $bam(glob("$bam_path/*.bam")){
    my $bamInfo = bamInfo($bam,$path_cleaned_reads,$total_ref_bases,$vcf_path);
  }

  return 0;
}

sub bamInfo{
  my($bam,$path_cleaned_reads,$total_ref_bases,$vcf_path)=@_;
  # Example filename: bam/sample1.fastq.gz-reference.sorted.bam
  # Example basename: sample1
  my $basename = basename($bam,qw(.fastq.gz-reference.sorted.bam));
			   
  print ("\n------------------------------------------------------------\n");
  print "\n\n\n---> Processing file: $basename...\n";
  print ("\n\n===================== QC on BAM File ===============\n");

  my $total_mapped_reads = `samtools view -c -F 4 $bam`; #Fetching total number of mapped reads  #NOTE: saving system() output to perl variable, returns 0; it sets the exit status to perl variable; Hence use ``
  my $mapped_reads = `samtools view -q 30 $bam | wc -l`;  #Fetching number of mapped reads with quality higher than 30
  my $mapped_reads_lt_30 = $total_mapped_reads - $mapped_reads; #Number of reads < Q30
               
  #Calculating total number of reads from clean_reads files
  my $file = "$path_cleaned_reads"."/"."$basename*gz";          
  my $line_count= `zcat $file | wc -l`;   #gettting total number of lines in fastq files
  my $tot_reads = $line_count/4;          #Dividing total number of lines by 4, to get total number of reads
  my $not_mapped = $tot_reads - $total_mapped_reads; # Calculating number of reads that are not mapped
  my $perc_mapped_reads = ($total_mapped_reads/$tot_reads) * 100; #Calculating percentage of mapped reads
  my $perc_mapped_reads_gt_30 = ($mapped_reads/$tot_reads) * 100; #Calculating percentage of mapped reads > Q30
			 
  print "\nTotal number of reads in File -> $basename: $tot_reads\n";  # Printing total number of reads
  print "Total number of mapped reads: $total_mapped_reads";
  print "Number of Reads Mapped > Q30: $mapped_reads";
  print "Number of Reads mapped < Q30: $mapped_reads_lt_30\n";
  print "Number of Reads not mapped: $not_mapped\n";
  print "Percent mapped reads: $perc_mapped_reads\n";
  print "Percent mapped reads > Q30: $perc_mapped_reads_gt_30\n";

  my $str1 = `samtools stats $bam | head -37 | grep 'reads properly paired:'`;
  $str1 =~s/^SN\s+//;        #Formatting output
  $str1 =~s/#.*//;
  print $str1;                         #reads properly paired:  ###### is stored as one string. Inorder to get the number, Splitting at ':'

  my @temp = split ":", $str1; 
  my $proper_paired = $temp[1];        #second element of the array is the number
  $proper_paired =~s/\s+//;           #Removing the white spaces before the number
  my $perc_proper_paired = ($proper_paired/$mapped_reads) * 100;



  my $str2 = `samtools stats $bam | head -37 | grep 'inward oriented pairs:'`;
  $str2 =~s/^SN\s+//;                        #Formatting output
  print $str2;

  my $str3 = `samtools stats $bam | head -37 | grep 'outward oriented pairs:'`;
  $str3 =~s/^SN\s+//;                        #Formatting output
  print $str3;

  my $str4 = `samtools stats $bam | head -37 | grep 'pairs with other orientation:'`;
  $str4 =~s/^SN\s+//;                        #Formatting output
  print $str4;


  print "Percentage of properly paired reads = $perc_proper_paired\n";


  my $str5 = `samtools stats $bam | head -37 | grep 'insert size average:'`;
  $str5 =~s/^SN\s+//;                        #Formatting output
  print $str5;


  my $baseslt20 = `samtools depth $bam | awk '( (\$3 >= 0) && (\$3 <= 19) )'| wc -l`;
  print "Number of bases with mapping coverage less than 20: $baseslt20";


  my $perc_reference_coverage = ((($total_ref_bases - $baseslt20)/$total_ref_bases)*100);
  print "Percent reference coverage = $perc_reference_coverage\n";



  my $tot_coverage_bases = `samtools depth $bam | cut -f 3| awk '{sum += \$_} END {print sum}'`;  # Total coverage of all bases
  print "Total coverage of bases (based off 20x at each position) = $tot_coverage_bases";

  my $average_coverage = $tot_coverage_bases/$total_ref_bases;
  print "Average Coverage = $average_coverage";


  print "\n\n\n=============================== QC on VCF File =========================\n";
  my $vcf_filename= "$vcf_path"."/".$basename."*.gz";


  print("Filters applied for sites that do not pass -
  a. DP20 = Depth is less than 20, the user-set coverage threshold
  b. RF0.95 = Reference variant consensus is less than 0.95, the user-set threshold
  c. AF 0.95 = Allele variant consensus is less than 0.95, the user-set threshold
  d. isIndel = Indels are not used for analysis in Lyve-SET
  e. masked = This site was masked using a bed file or other means
  f. str10 = Less than 10% or more than 90% of variant supporting reads on one strand
  g. indelError = Likey artifcat due to indel reads at this position\n\n");
               
  my $filter_codes=`zcat $vcf_filename | grep -v '^#' | cut -f 7 | sort | uniq -c`;
  print("$filter_codes");

}


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

  #Declare a new array for the sorting the unsorted array. Use $$ to reference the other array. 
  #Sort the array within an array on $num (the second element in the unsorted array) in ascending order (a <=> b) on GC content 
  #my @sorted= sort {$$a[2] <=> $$b[2]} @unsorted;
               
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
  while (my $seq = $seq_in->next_seq ()) 
                { 

  #Initialize variables; start the counter at 0 for each contig that comes through the loop  

  $id= $seq->id; #Get the header Id of each contig
  $length +=$seq->length;  #Declare variable for finding the length of each contig sequence

  $string = $seq->seq;   #Each contig seq is treated as a string


  }

  return $length;

}

sub usage{
  "Runs a diagnosis on a finished Lyve-SET run
  Usage: ./hqSNP_diagnose.pl projectdir
  "
}
