#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse dirname basename/;
use File::Temp qw/tempdir tmpfile/;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help outdir=s format=s));
  die usage() if($$settings{help});
  $$settings{format}||="tsv"; # by default, input format is tsv

  # Get the output directory and spreadsheet, and make sure they exist
  $$settings{outdir}||=die "ERROR: need outdir parameter\n".usage();
  mkdir $$settings{outdir} if(!-d $$settings{outdir});
  my $tsv=$ARGV[0] || die "ERROR: need tsv file!\n".usage();
  die "ERROR: cannot find $tsv" if(!-e $tsv);

  # Read the spreadsheet
  my $infoTsv=readTsv($tsv,$settings);

  # download everything
  downloadEverything($infoTsv,$settings);

  return 0;
}

sub readTsv{
  my($tsv,$settings)=@_;

  # For the fastq-dump command
  my $seqIdTemplate='@$ac_$sn[_$rn]/$ri';
  
  my $d={}; # download hash
  my $have_reached_biosample=0; # marked true when it starts reading entries
  my @header=(); # defined when we get to the biosample_acc header row
  open(TSV,$tsv) or die "ERROR: could not open $tsv: $!";
  while(<TSV>){
    s/^\s+|\s+$//g; # trim whitespace
    next if(/^$/);

    ## read the contents
    # Read biosample rows
    if($have_reached_biosample){
      my $tmpdir=tempdir("$0XXXXXX",TMPDIR=>1,CLEANUP=>1);
      # Get an index of each column
      my %F;
      @F{@header}=split(/\t/,$_);

      # SRA download command
      if($F{srarun_acc}){
        $$d{$F{srarun_acc}}{download}="fastq-dump --defline-seq '$seqIdTemplate' --defline-qual '+' --split-files -O $tmpdir --gzip $F{srarun_acc} ";
        $$d{$F{srarun_acc}}{name}=$F{strain} || die "ERROR: $F{srarun_acc} does not have a strain name!";
        $$d{$F{srarun_acc}}{type}="sra";
        $$d{$F{srarun_acc}}{tempdir}=$tmpdir;
      }

      # GenBank download command
      if($F{genbankassembly}){
        $$d{$F{genbankassembly}}{download}="efetch -format gbwithparts -db nuccore -id $F{genbankassembly} > $tmpdir/$F{genbankassembly}.gbk";
        $$d{$F{genbankassembly}}{name}=$F{strain} || die "ERROR: $F{genbankassembly} does not have a strain name!";
        $$d{$F{genbankassembly}}{type}="genbank";
        $$d{$F{genbankassembly}}{tempdir}=$tmpdir;
      }

    } elsif(/^biosample_acc/){
      $have_reached_biosample=1;
      @header=split(/\t/,lc($_));
      next;
    }
    # metadata
    else {
      my ($key,$value)=split /\t/;
      $key=lc($key);
      $$d{$key}=$value;
    }

  }
  close TSV;

  #print Dumper $d;die;
  return $d;
}

sub downloadEverything{
  my($d,$settings)=@_;

  # Read each entry one at a time.  Each entry is a hash
  # consisting of: type, name, download, tempdir.
  while(my($key,$value)=each(%$d)){
    # Only download entries which are hash values
    next if(ref($value) ne "HASH" || !defined($$value{download}));

    # Get some local variables to make it more readable downstream
    my($type,$name,$download,$tempdir)=($$value{type},$$value{name},$$value{download},$$value{tempdir});

    # Run the download command found in the entry.
    # The actual download command is found in the if statements
    # in case we need to debug and focus on a specific statement.
    logmsg "Downloading $name/$type to $tempdir";

    # Rename/move the downloaded files appropriately.
    if($type eq 'sra'){
      # Perform the download
      system($download);
      die "ERROR downloading $name/$type" if $?;
      # Move fastq.gz files to $strain.$SRR_1.fastq.gz.
      for my $file(glob("$tempdir/*.fastq.gz")){
        my $b=basename($file);
        system("mv -v $file $$settings{outdir}/$name.$b");
        die "ERROR moving $file to $name.$b" if $?;
      }
    } elsif($type eq 'genbank'){
      # Perform the download
      system($download);
      die "ERROR downloading $name/$type" if $?;
      # Move gbk files to $strain.gbk.
      # Only one gbk should be present. Others will be lost,
      # if they exist.
      for my $file(glob("$tempdir/*.gbk")){
        my $b=basename($file);
        system("mv -v $file $$settings{outdir}/$name.gbk");
        die "ERROR moving $file to $name.gbk" if $?;
      }
    } else {
      die "INTERNAL ERROR: do not know how to deal with type $type";
    }
  }

  # I can't think of any useful return value at this time.
  return 1;
}

sub usage{
  "Reads a standard dataset spreadsheet and downloads its data
  Usage: $0 -o outdir spreadsheet.dataset.tsv
  --format     tsv    The input format. Default: tsv. No other format
                      is accepted at this time.
  "
}

