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
  GetOptions($settings,qw(help outdir=s));
  die usage() if($$settings{help});

  # Get the output directory and spreadsheet, and make sure they exist
  $$settings{outdir}||=die "ERROR: need outdir parameter\n".usage();
  mkdir $$settings{outdir} if(!-d $$settings{outdir});
  my $tsv=$ARGV[0] || die "ERROR: need tsv file!\n".usage();
  die "ERROR: cannot find $tsv" if(!-e $tsv);

  # Read the spreadsheet
  my $info=readSpreadsheet($tsv,$settings);

  # download everything
  downloadEverything($info,$settings);

  return 0;
}

sub readSpreadsheet{
  my($tsv,$settings)=@_;

  # For the fastq-dump command
  my $seqIdTemplate='@$sn[_$rn]/$ri';
  
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
      $$d{$F{srarun_acc}}{download}="fastq-dump --defline-seq '$seqIdTemplate' --defline-qual '+' --split-files -O $tmpdir --gzip $F{srarun_acc} ";
      $$d{$F{srarun_acc}}{name}=$F{strain};
      $$d{$F{srarun_acc}}{type}="sra";

      # GenBank download command
      #(for i in `seq 1 7`; do id="JASV0100000$i"; efetch -db nucleotide -format fasta -id $id; done;) | grep . | sha256sum –
      $$d{$F{genbankassembly}}{download}="efetch -format gbwithparts -db nuccore -id $F{genbankassembly} > $tmpdir/$F{genbankassembly}.gbk";
      $$d{$F{genbankassembly}}{name}=$F{strain};
      $$d{$F{genbankassembly}}{type}="genbank";

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
  return $d;
}

sub downloadEverything{
  my($d,$settings)=@_;

  while(my($key,$value)=each(%$d)){
    next if(!defined($$value{download}));

    my($type,$name,$download)=($$value{type},$$value{name},$$value{download});

    logmsg "Downloading $name/$type";
    system($download);
    die "ERROR downloading $name/$type" if $?;

    # Download things appropriately and rename them.
    if($type eq 'sra'){
      for my $file(glob("*.fastq.gz")){
        my $b=basename($file);
        system("mv $file $$settings{outdir}/$name.$b");
        die "ERROR moving $file to $name.$b" if $?;
      }
    } elsif($type eq 'genbank'){
      for my $file(glob("*.gbk")){
        my $b=basename($file);
        system("mv $file $$settings{outdir}/$name.gbk");
        die "ERROR moving $file to $name.gbk" if $?;
      }
    }
  }
}

sub usage{
  "Reads a standard dataset spreadsheet and downloads its data
  Usage: $0 -o outdir spreadsheet.dataset.tsv
  "
}

