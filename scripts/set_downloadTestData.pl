#!/usr/bin/env perl
# Download all the test data for Lyve-SET
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Cwd qw/realpath/;
use File::Basename qw/basename/;
use File::Copy qw/copy move/;
use File::Slurp qw/read_file/;
use File::Temp qw/ tempfile tempdir /;
use Bio::Perl;

use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH}="$ENV{PATH}:$FindBin::RealBin/../lib/edirect";

# get project test data
my $testDir="$FindBin::RealBin/../testdata";
my @test=_uniq(map({basename($_) if(-d $_)} glob("$testDir/*")));
my $testdata=join(", ",@test);

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{tempdir}||=tempdir(CLEANUP=>1);

  # put all the datasets in lowercase
  my @dataset=@ARGV;
  $_=lc($_) for(@dataset);
  # index all the datasets
  my %dataset;
  @dataset{@dataset}=@dataset;
  # if 'all' is invoked, then capture all the databases
  if($dataset{all}){
    @dataset=@test;
    @dataset{@dataset}=@dataset;
  }

  # Download the datasets
  for(@dataset){
    logmsg "Downloading $_";
    downloadDataset($_,$settings);
  }

  return 0;
}

sub downloadDataset{
  my($dataset,$settings)=@_;
  
  my $testdir="$testDir/$dataset";

  # check for files
  my @readsLocal=glob("$testdir/reads/*.fastq.gz");
  my @assemblyLocal=glob("$testdir/asm/*.fasta");
  # assembly files?
  my @readsRemote=read_file("$testdir/SRA",{err_mode=>"quiet"});
  my @assemblyRemote=read_file("$testdir/NUCLEOTIDE",{err_mode=>"quiet"});

  my @reads=(@readsLocal,@readsRemote);
  my @assembly=_uniq(@assemblyLocal,@assemblyRemote);

  die "ERROR: no reads found in $testdir\n" if(!@reads);
  die "ERROR: no reference genome found found in $testdir/NUCLEOTIDE" if(!@assembly);
  chomp(@reads,@assembly,@readsRemote,@assemblyRemote);

  # get the reads
  mkdir("$testdir/reads");
  for(@readsRemote){
    $_=~s/^\s+|\s+$//g; # trim
    $_=~s/#.*//;        # ignore comments
    next if(/^$/);      # if there's nothing left on this line anymore, then skip it

    my ($r,$realname)=split(/\t/,$_);
    $realname=$r if(!$realname);
    $_=~s/^\s+|\s+$//g for($r,$realname);
    downloadReads($r,$realname,$testdir,$settings);
  }

  # assemblies: the first assembly is the reference genome
  mkdir("$testdir/asm");
  for(@assemblyRemote){
    $_=~s/^\s+|\s+$//g; # trim
    $_=~s/#.*//;        # ignore comments
    next if(/^$/);      # if there's nothing left on this line anymore, then skip it

    my ($r,$realname)=split(/\t/,$_);
    $realname=$r if(!$realname);
    $_=~s/^\s+|\s+$//g for($r,$realname);
    _downloadAssembly($r,$realname,$testdir,$settings);
  }

  # reference genome
  mkdir("$testdir/reference");
  my ($r,$realname)=split(/\t/,$assemblyRemote[0]);
  $realname=$r if(!$realname);
  $_=~s/^\s+|\s+$//g for($r,$realname);
  logmsg "Establishing the first nucleic acid sequence $r ($realname) as the reference genome";
  copy("$testdir/asm/$realname.fasta","$testdir/reference");
  die "ERROR: could not copy $realname.fasta to $testdir/reference: $!" if $?;

  return 1;
}

sub downloadReads{
  my($SRR,$realname,$testdir,$settings)=@_;
  my $finalFile="$testdir/reads/$realname.fastq.gz";
  if(-e $finalFile){
    logmsg "$finalFile found. Not downloading again";
    return $finalFile;
  }

  # Download using SRA toolkit. 
  # Prefetch knows whether a file has been downloaded already and will skip it as warranted
  # So there is no need for a conditional statement here
  logmsg "Downloading $SRR ($realname) with prefetch command";
  system("prefetch $SRR");
  die "ERROR with the prefetch command. Is the sra-toolkit installed?  Command was\n  prefetch $SRR" if $?;

  # uncompress
  my $fastqDumpDir="$$settings{tempdir}/$SRR.fastq";
  logmsg "Converting sra file to fastq with fastq-dump";
  if(! -d $fastqDumpDir){
    my $deflineArg='@$sn[_$rn]/$ri';
    my $command="fastq-dump -v -v -v -v -v --defline-seq '$deflineArg' --defline-qual '+' --split-files -O $fastqDumpDir.tmp $SRR && mv $fastqDumpDir.tmp $fastqDumpDir";
    system($command);
    die "ERROR with fastq-dump\n  $command" if $?;
  }

  # shuffle reads
  logmsg "Shuffling the fastq reads into the Lyve-SET directory";
  if(!-e $finalFile){
    system("run_assembly_shuffleReads.pl $fastqDumpDir/${SRR}_1.fastq $fastqDumpDir/${SRR}_2.fastq | gzip -c > $finalFile.tmp && mv -v $finalFile.tmp $finalFile");
    die "ERROR with either run_assembly_shuffleReads.pl or gzip" if $?;
    move("$finalFile.tmp",$finalFile);
    die "ERROR renaming $finalFile.tmp to $finalFile: $!" if $?;
  }

  system("rm -rvf $fastqDumpDir");
  logmsg "WARNING: error in removing temporary files" if $?;

  return $finalFile;
}

sub _downloadAssembly{
  my($NC,$realname,$datadir,$settings)=@_;
  my $tmpFile="$$settings{tempdir}/$NC.fasta";
  my $filteredFile="$datadir/asm/$realname.fasta";
  return $filteredFile if(-f $filteredFile);

  my $cmd="esearch -db nucleotide -query $NC | efetch -format fasta > $tmpFile";
  logmsg "Downloading $NC. Command:\n  $cmd";
  system($cmd);
  die "ERROR with efetch or esearch" if $?;
  die "ERROR: downloaded a zero-byte file for $NC. Is it the correct identifier?\n" if (-s $tmpFile < 1);

  # process the downloaded fasta file
  my $numContigsPassed=0;
  my $in=Bio::SeqIO->new(-format=>"fasta",-file=>$tmpFile);
  my $out=Bio::SeqIO->new(-format=>"fasta",-file=>">$filteredFile.tmp");
  while(my $seq=$in->next_seq){
    # filter the sequence to contigs >500bp
    next if($seq->length < 500);
    # make the defline kosher
    my $id=$seq->id;
    $id=~s/\:/_/g;
    $seq->id($id);
    $seq->desc("descriptionRemovedByLyveSet");

    # Write it all out
    $out->write_seq($seq);
    $numContigsPassed++;
  }
  $in->close;
  $out->close;
  logmsg "$numContigsPassed contigs passed the filter for $NC";
  die "ERROR: no contigs passed the filter!" if($numContigsPassed < 1);

  move("$filteredFile.tmp",$filteredFile);
  die "Error with moving $filteredFile.tmp to $filteredFile" if $?;

  return $filteredFile;
}

# Return unique and non-empty array of elements
sub _uniq{
  my %seen;
  my @arr=grep { !$seen{$_}++} @_;
  my @nonEmpty;
  for(@arr){
    push(@nonEmpty,$_) if($_);
  }
  return @nonEmpty;
}

sub usage{
  "Download all the test data under the testdata directory
  Usage: $0 data_set
  data_set can be either 'all' or one of the following:
    $testdata
  All data are downloaded to $testDir
  "
}
