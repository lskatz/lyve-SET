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
  GetOptions($settings,qw(help do-nothing numcpus=i fast)) or die $!;
  $$settings{numcpus}||=1;

  my ($dataset,$project,@setArgv)=@ARGV;
  die usage() if(!@ARGV || $$settings{help});
  die "ERROR: need a dataset name\n".usage() if(!$dataset);
  push(@setArgv,"--fast") if($$settings{fast});
  $project||=$dataset;
  #push(@setArgv,qw(--rename-taxa s/\\\\\\..*//)); # doesn't seem to put in the backslashes into the final processMsa command yet

  if(!-d $project){
    getData($dataset,$project,$settings);
  }
  launchSet($project,\@setArgv,$settings);

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

  # link the reference genome assembly
  for(glob("$datadir/$dataset/asm/*.fasta")){
    command("set_manage.pl $project --change-reference $_",$settings);
  }

  # All steps after this are only gravy and so skip them
  # if we are trying to be fast.
  return if($$settings{fast});

  #link the assemblies
  for(glob("$datadir/$dataset/asm/*.fasta")){
    command("set_manage.pl $project --add-assembly $_",$settings);
  }
}

sub launchSet{
  my($project,$setArgv,$settings)=@_;
  my $setArgvStr=join(" ",@$setArgv);
  command("launch_set.pl $project --numcpus $$settings{numcpus} $setArgvStr",$settings);
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
  Usage: $0 dataset project
  dataset names could be one of the following:\n    ".join(", ",@dataname)."
  NOTE: project is the output directory for Lyve-SET

  --numcpus 1  # How many cpus you want to use
  --do-nothing # To print the commands but do not run system calls
  --fast       # Use some shortcuts: do not simulate reads; use --fast for launch_set.pl
  -- [......]  # Put any parameters for launch_set.pl after a double dash and a space
               # If used, then the user must specify [project]
  "
}
