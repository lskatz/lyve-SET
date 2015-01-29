#!/usr/bin/env perl
# Tests SET
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename/;
use FindBin;

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}

# find out what the test data could be
my $datadir="$FindBin::RealBin/../testdata";
my (@dataname,%dataname);
for my $path(glob("$datadir/*")){
  my $name=basename($path);
  push(@dataname,$name) if($name=~/^.+$/ && -d $path);
  $dataname{$name}=1;
}

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help do-nothing numcpus=i)) or die $!;
  $$settings{numcpus}||=1;

  my ($dataset,$project)=@ARGV;
  die "ERROR: need a dataset name\n".usage() if(!$dataset);
  $project||=$dataset;

  getData($dataset,$project,$settings);
  launchSet($project,$settings);

  return 0;
}

sub getData{
  my($dataset,$project,$settings)=@_;

  # download the data
  command("set_downloadTestData.pl $dataset",$settings);

  # create the project
  command("set_manage.pl --create $project",$settings);

  # link the reads
  for(glob("$datadir/$dataset/reads/*.fastq.gz")){
    command("set_manage.pl $project --add-reads $_",$settings);
  }

  #link the assemblies
  for(glob("$datadir/$dataset/asm/*.fasta")){
    command("set_manage.pl $project --add-assembly $_",$settings);
  }

  # link the reference genome assembly
  for(glob("$datadir/$dataset/asm/*.fasta")){
    command("set_manage.pl $project --change-reference $_",$settings);
  }
}

sub launchSet{
  my($project,$settings)=@_;
  command("launch_set.pl $project --numcpus $$settings{numcpus} --numnodes 50",$settings);
}

sub command{
  my($cmd,$settings)=@_;
  logmsg "COMMAND\n  $cmd";
  my $exit_code=system($cmd) if(!$$settings{'do-nothing'});
  if($exit_code){
    logmsg "ERROR: died while trying to run COMMAND\n  $cmd";
    logmsg "See if you have the correct software installed for this command.";
    die $exit_code;
  }
  return $exit_code;
}

sub usage{
  "Runs a test dataset with Lyve-SET
  Usage: $0 dataset [project]
  dataset names could be one of the following:\n    ".join(", ",@dataname)."
  NOTE: project will be the name of the dataset, if it is not given

  --numcpus 1  How many cpus you want to use
  --do-nothing To print the commands but do not run system calls
  "
}
