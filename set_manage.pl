#!/usr/bin/env perl
# Manage a Lyve-SET project
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Cwd qw/realpath/;
use File::Basename qw/basename/;
use File::Copy qw/copy/;

sub logmsg{print STDERR "@_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help create delete add-reads=s remove-reads=s add-assembly=s remove-assembly=s change-reference=s)) or die $!;
  die usage() if($$settings{help});
  die "ERROR: need a SET project\n".usage() if(!@ARGV);
  my $project=shift(@ARGV);

  createProjectDir($project,$settings) if($$settings{create});
  if(!is_project($project,$settings)){
    die "$project is not a SET project!";
  }
  addReads($project,$settings) if($$settings{'add-reads'});
  addAssembly($project,$settings) if($$settings{'add-assembly'});
  changeReference($project,$settings) if($$settings{'change-reference'});
  removeAssembly($project,$settings) if($$settings{'remove-assembly'});
  removeReads($project,$settings) if($$settings{'remove-reads'});
  deleteProject($project,$settings) if($$settings{delete});

  return 0;
}

sub createProjectDir{
  my($dir,$settings)=@_;
  die "ERROR: directory $dir already exists!" if(-d $dir);
  mkdir $dir;
  die $! if $?;
  for(qw(vcf vcf/unfiltered msa bam reads reference tmp asm)){
    mkdir "$dir/$_";
    die $! if $?;
  }
  return $dir;
}

sub is_project{
  my($project,$settings)=@_;
  for(qw(vcf vcf/unfiltered msa bam reads reference tmp asm)){
    return 0 if(!-e "$project/$_");
  }
  return 1;
}

sub addReads{
  my($project,$settings)=@_;
  my $reads=File::Spec->rel2abs($$settings{'add-reads'});
  my $symlink="$project/reads/".basename($reads);
  symlink($reads,$symlink);
  logmsg "$reads => $symlink";
  return 1;
}
sub addAssembly{
  my($project,$settings)=@_;
  my $asm=File::Spec->rel2abs($$settings{'add-assembly'});
  my $symlink="$project/asm/".basename($asm);
  symlink($asm,$symlink);
  logmsg "$asm => $symlink";
  return 1;

}
sub removeReads{
  my($project,$settings)=@_;
  my $name=$$settings{'remove-reads'};
  for my $file("$project/reads/$name",glob("$project/bam/$name*"),glob("project/vcf/$name*"),glob("project/vcf/unfiltered/$name*")){
    unlink($file) or logmsg "Warning: could not remove $file";
  }
  return 1;
}
sub removeAssembly{
  my($project,$settings)=@_;
  my $name=$$settings{'remove-assembly'};
  for my $file("$project/asm/$name"){
    unlink($file) or logmsg "Warning: could not remove $file";
  }
  return 1;
}

sub changeReference{
  my($project,$settings)=@_;
  my $ref=$$settings{'change-reference'};
  my $newPath="$project/reference/".basename($ref);
  copy($ref,$newPath) or die "ERROR: could not copy $ref to $newPath";
  logmsg "cp $ref => $newPath";

  unlink("$project/reference/reference.fasta");
  symlink(basename($newPath),"$project/reference/reference.fasta");
  logmsg "symlink $newPath => $project/reference/reference.fasta";
  return 1;
}

sub deleteProject{
  my($project,$settings)=@_;
  system("rm -rfv '$project'");
  die if $?;
  return 1;
}

sub usage{
  "Manages a Lyve-SET project and confirms that a directory is a Lyve-SET project.
  Usage: $0 setProject [options]
  --create Create a set project
  --delete Delete a set project
  --add-reads file.fastq.gz       Add reads to your project (using symlink)
  --remove-reads file.fastq.gz    Remove reads from your project (using rm)
  --add-assembly file.fasta       Add an assembly to your project (using symlink)
  --remove-assembly file.fasta    Remove an assembly from your project (using rm)
  --change-reference file.fasta   Add a reference assembly to your project and remove any existing one (using cp)
  Examples
    set_manage.pl projDir --add-reads ../../reads/unknown/file.cleaned.fastq.gz
    set_manage.pl projDir --remove-reads file.cleaned.fastq.gz
  "
}
