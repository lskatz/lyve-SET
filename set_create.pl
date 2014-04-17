#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  die usage() if($$settings{help});
  my @dir=@ARGV;
  die "ERROR: need project directory(ies)\n".usage() if(!@dir);

  for my $dir(@dir){
    createProjectDir($dir,$settings);
    logmsg "Created project directory $dir";
  }
  return 0;
}

sub createProjectDir{
  my($dir,$settings)=@_;
  die "ERROR: project $dir already exists!" if(-d $dir);
  mkdir $dir;
  die $! if $?;
  for(qw(vcf vcf/unfiltered msa bam reads reference tmp)){
    mkdir "$dir/$_";
    die $! if $?;
  }
  return $dir;
}

sub usage{
  "Creates a directory structure suitable for SET
  usage: $0 projectDir
  "
}
