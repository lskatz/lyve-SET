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
  my $settings={importasm=>1};
  GetOptions($settings,qw(help do-nothing numcpus=i numnodes=i importasm! fast)) or die "$!\n  Please remember, if you are trying to use an option for launch_set.pl, use the '--' flag first. See the usage statement for more details.";
  $$settings{numcpus}||=1;
  $$settings{numnodes}||=50;

  my ($dataset,$project,@setArgv)=@ARGV;
  die usage() if(!@ARGV || $$settings{help});
  die "ERROR: need a dataset name\n".usage() if(!$dataset);
  $project||=$dataset;

  # Put all 'fast' options here
  if($$settings{fast}){
    $$settings{numnodes}=9999 if($$settings{numnodes}<9999);
    $$settings{importasm}=0;
    push(@setArgv,"--fast");
  }
  
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
  return if(!$$settings{importasm});

  #link the assemblies
  for(glob("$datadir/$dataset/asm/*.fasta")){
    command("set_manage.pl $project --add-assembly $_",$settings);
  }
}

sub launchSet{
  my($project,$setArgv,$settings)=@_;
  my $setArgvStr=join(" ",@$setArgv);
  command("launch_set.pl $project --numcpus $$settings{numcpus} --numnodes $$settings{numnodes} $setArgvStr",$settings);
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
  Usage: $0 dataset project [-- --launch_set.pl options ...]
  dataset names could be one of the following:\n    ".join(", ",@dataname)."
  NOTE: project is the output directory for Lyve-SET

  --numcpus 1   # How many cpus you want to use
  --numnodes 20 # How many maximum nodes to use
  --do-nothing  # To print the commands but do not run system calls
  --noimportasm # Do not import assmblies and simulate reads
  --fast        # Same as: --noimportasm --numnodes 9999 -- --fast
  -- [......]   # Put any parameters for launch_set.pl after a double dash and a space
  "
}
