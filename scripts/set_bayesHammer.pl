#!/usr/bin/env perl

## Simple version of the BayesHammer wrapper script.
## Assumes user provides only paired-end Illumina reads as input.
## Paired reads are provided as output, Unpaired reads discarded.
## Requires SPAdes 3.0 or higher and cg_pipeline/scripts/run_assembly_shuffleReads.pl


use strict;
use warnings;
use Data::Dumper;
use Cwd;
use Getopt::Long;
use Bio::Perl;
use File::Basename qw/basename fileparse/;
use File::Temp qw/tempdir/;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

use LyveSET qw/logmsg/;

local $0=basename($0);
exit (main());

sub main{
  my $settings={};
  
  die usage() if(@ARGV < 2);  
  GetOptions($settings,qw(help tmpdir=s numcpus=s)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("XXXXXX",TMPDIR=>1,CLEANUP=>1);
  die usage() if($$settings{help});
  logmsg "Temporary directory is $$settings{tempdir}";

  my @infile=@ARGV;
  for my $infile(@infile){
    errorCorrect($infile,$settings);

  }
  
  return 0;
}

sub errorCorrect{
  my($infile,$settings)=@_;

  # Is it paired end?
  my $spadesInParam="--12 $infile";
  if(`run_assembly_isFastqPE.pl` == 0){
    $spadesInParam="-s $infile";
  }

  # Run spades
  system("spades.py $spadesInParam -o $$settings{tempdir} --only-error-correction --disable-gzip-output -t $$settings{numcpus} >&2");
  die "ERROR with spades.py" if $?;

  
}
   my $R1_head = IsFastq($reads1);
   my $R2_head = IsFastq($reads2);  
 
  ## extract R1 strand identifier only
   my @Fastq_head = split(' ', $R1_head);
   $R1_head = $Fastq_head[1];

  ## extract R2 strand identifier only
   @Fastq_head = split(' ', $R2_head);
   $R2_head = $Fastq_head[1];

  # print $R1_head, " ", $R2_head, "\n";  

   my $Folder = '';

   if($$settings{tmpdir})
     {
        $Folder = BuildSubFolder($$settings{tmpdir});
     }
   else
    {
        $Folder = BuildSubFolder("Temp");
    }
   
   RunBayesHammer($reads1, $reads2, $Folder, $threads);
  
   ReadBayesFolder($Folder, $R1_head, $R2_head);

   ## search *.cor.fastq files in $Folder and shuffle them
   OutputShuffledReads($Folder);

   return 0;
}

sub IsFastq{

    my($gzinfile) = @_;
    
    my $first_ID = '';

    if($gzinfile =~ /\.fastq\.gz\b/)
      {
         logmsg("reading $gzinfile to obtain strand identifier");
 	 my $infile = substr($gzinfile, 0, length($gzinfile) - 3);
         gunzip $gzinfile => $infile or die "gunzip failed: $GunzipError\n";
         open IN, $infile || die "Can't find $infile $!";
         $first_ID = <IN>;
         logmsg("removing unzipped $infile");
         system("rm $infile");
        
      }
    elsif($gzinfile =~ /\.fastq\b/)
    {
      open IN, $gzinfile || die "Can't find $gzinfile $!";
      $first_ID = <IN>;
      
    }
    else
     {
	 die "File, $gzinfile, does not have valid fastq format.\n $!";
         #.usage();
     }

        return $first_ID;
}

sub BuildSubFolder{

    my($infile) = @_;
   my $SeqName = $infile;
   $SeqName =~ s/\.fastq(\.gz)?//g;
   $SeqName =~ s/_R(1|2)(_001)?//g; 
   $SeqName = "BayesHam_".$SeqName;
   die "ERROR: directory $SeqName already exists." if (-d $SeqName); 
   mkdir $SeqName;
   logmsg("mkdir $SeqName");
   return $SeqName;
}

sub RunBayesHammer{
    
    my($infile1, $infile2, $Folder, $threads) = @_;
   
    logmsg("Starting SPAdes error correction\n"); 
    system("spades.py -1 $infile1 -2 $infile2 -o $Folder --only-error-correction --disable-gzip-output -t $threads > $Folder/STDERR");
    
    return 0;
}

sub ReadBayesFolder{
    
    my($Folder, $R1, $R2) = @_;
   
    my $a = cwd();
   
    opendir(DIR, $Folder."/corrected") or die "Folder, $Folder/corrected, does not exist$!"; 

   while (my $file = readdir(DIR))
    {
	 #print $file, "\n";
       
        if($file =~ /_R1_.*fastq\b/)
        {
	  logmsg( "$file, in R1!" );
          open(FWD, $Folder."/corrected/".$file) || die "File, $file, cannot be opened $!";
          my $ofile = $file;
             $ofile =~ s/\.00\.0_0//g;
          open OUTFWD, ">$Folder/$ofile" or die $!;
          
          while(<FWD>)
           {
             if(($. % 4) == 1)
               {
                 $_ =~ s/\n/ /g;
                 print OUTFWD $_, $R1, "\n";   
               }
             elsif(($. % 4) == 3)
               {
                 print OUTFWD substr($_, 0, 1), "\n";  
               } 
             else
              {
                print OUTFWD $_;
              }
          }
	}
      elsif($file =~ /_R2_.*fastq\b/)
        {
	   logmsg( "$file, in R2!\n" );
          open(REV, $Folder."/corrected/".$file) || die "File, $file, cannot be opened $!";
          my $ofile = $file;
	   $ofile =~ s/\.00\.0_0//g;
          open OUTREV, ">$Folder/$ofile" or die $!;
          
          while(<REV>)
           {
             if(($. % 4) == 1)
               {
                $_ =~ s/\n/ /g;
                print OUTREV $_, $R2, "\n";   
               }
             elsif(($. % 4) == 3)
               {
                 print OUTREV substr($_, 0, 1), "\n";  
               } 
             else
              {
         	 print OUTREV $_;
              }
           } 
       }    
     }    

  ## Clean up subfolders ##
   
    logmsg("removing files in $Folder/corrected");
    system("rm -r $Folder/corrected/*");

    logmsg("removing $Folder/input_dataset.yaml");
    system("rm $Folder/input_dataset.yaml");
  
    logmsg("removing $Folder/corrected and $Folder/tmp");
    system("rmdir $Folder/corrected");
    system("rmdir $Folder/tmp");

    closedir(DIR);

    return 0;
}


sub OutputShuffledReads{
    
    my($Folder) = @_;
   
    my $a = cwd();
   
    opendir(DIR, $Folder) or die "Folder, $Folder, does not exist$!"; 
 
    my $Fastq_pair = '';
    my $tmp_file_str = '';    
    while (my $file = readdir(DIR))
     {
      if($file =~ /.*R1.*fastq\.cor\.fastq/)
        {
	    $Fastq_pair = $Folder."/".$file." ";
        }
      elsif($file =~ /.*R2.*fastq\.cor\.fastq/)
        {
           $tmp_file_str = $Folder."/".$file." ";
        }   
     } 

     $Fastq_pair = $Fastq_pair.$tmp_file_str;
    
    ## output shuffled reads to STDIN
    logmsg("run_assembly_shuffleReads.pl $Fastq_pair");
    system("run_assembly_shuffleReads.pl $Fastq_pair");
    
    return 0;
}   

sub usage{
  "$0: Corrects Illumina reads using BayesHammer in SPAdes
  Usage: 
    $0 reads.shuffled.fastq[.gz] > cleaned.shuffled.fastq
    --numcpus 1
    --tempdir <optional>

  "
}
