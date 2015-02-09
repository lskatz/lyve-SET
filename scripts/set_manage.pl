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
use File::Copy qw/copy move/;
use File::Slurp qw/read_file/;
use Bio::Perl;

use FindBin;
use lib "$FindBin::RealBin/../lib";
$ENV{PATH}="$ENV{PATH}:$FindBin::RealBin/../lib/edirect";

# The directories a project should have
my @projectSubdir=qw(vcf vcf/unfiltered msa bam reads reference tmp asm log);

# get project test data
my $testDir="$FindBin::RealBin/../testdata";
my @test=_uniq(map({basename($_) if(-d $_)} glob("$testDir/*")));
my $testdata=join(", ",@test);

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help create test=s delete add-reads=s remove-reads=s add-assembly=s remove-assembly=s change-reference=s numnodes=i numcpus=i)) or die $!;
  die usage() if($$settings{help});
  die "ERROR: need a SET project\n".usage() if(!@ARGV);
  my $project=shift(@ARGV);
  $$settings{create}=1 if($$settings{test});

  # parameters for compatibility reasons
  $$settings{numnodes}||=1;
  $$settings{numcpus}||=1;

  createProjectDir($project,$settings) if($$settings{create});
  if(!is_project($project,$settings)){
    my(undef,$dirs)=is_project($project,$settings);
    die "$project is not a SET project! Missing $dirs";
  }
  $$settings{tempdir}||="$project/tmp";

  addSampleData($project,$settings) if($$settings{test});
  addReads($project,$settings) if($$settings{'add-reads'});
  addAssembly($project,$settings) if($$settings{'add-assembly'});
  changeReference($project,$settings) if($$settings{'change-reference'});
  removeAssembly($project,$settings) if($$settings{'remove-assembly'});
  removeReads($project,$settings) if($$settings{'remove-reads'});
  deleteProject($project,$settings) if($$settings{delete});
  removeGenome($project,$settings) if($$settings{'remove-genome'});

  return 0;
}

sub createProjectDir{
  my($dir,$settings)=@_;
  die "ERROR: directory $dir already exists!" if(-d $dir);
  mkdir $dir;
  die $! if $?;
  for(@projectSubdir){
    mkdir "$dir/$_";
    logmsg "mkdir $dir/$_";
    die $! if $?;
  }
  return $dir;
}

sub addSampleData{
  my($project,$settings)=@_;
  
  my $testdir="$testDir/$$settings{test}";

  # check for files
  my @reads=glob("$testdir/reads/*.fastq.gz");
  my @reference=glob("$testdir/reference/*.fasta");
  my @assembly=glob("$testdir/asm/*.fasta");

  die "ERROR: no reads found in $testdir\n" if(!@reads);
  die "ERROR: no reference found in $testdir" if(!@reference);
  logmsg "WARNING: no assemblies found in $testdir" if(!@assembly);
  chomp(@reads,@reference,@assembly);

  # get the reads
  for(@reads){
    local $$settings{'add-reads'}=$_;
    my $path=addReads($project,$settings);
  }

  # add assemblies
  for(@assembly){
    local $$settings{'add-assembly'}=$_;
    addAssembly($project,$settings);
  }

  # get the reference
  local $$settings{'change-reference'}=$reference[0];
  logmsg "Using ".$$settings{'change-reference'}." as the reference genome";
  changeReference($project,$settings);

  logmsg "Done! To test these data, please execute\n launch_set.pl $project";
  return 1;
}

sub is_project{
  my($project,$settings)=@_;
  my $dirs="";
  for(@projectSubdir){
    $dirs.="$project/$_ " if(!-e "$project/$_");
  }

  if($dirs){
    return (0,$dirs) if wantarray;
    return 0;
  }
  return (1,"") if wantarray;
  return 1;
}

sub addReads{
  my($project,$settings)=@_;
  my $reads=File::Spec->rel2abs($$settings{'add-reads'});
  my $newPath="$project/reads/".basename($reads);
  if(-e $reads){
    my $symlink=$newPath;
    symlink($reads,$symlink);
    logmsg "$reads => $symlink";
  } else {
    logmsg "NOTE: could not find file $reads";
    logmsg "I will see if it's on NCBI";
    my $b=basename($reads,qw(.sra .fastq .fastq.gz));
    $newPath="$project/reads/$b.fastq.gz";
    my $three=substr($b,0,3);
    my $six=substr($b,0,6);
    system("wget --continue -O $$settings{tempdir}/$b.sra 'ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$three/$six/$b/$b.sra'");
    die "ERROR with wget, or maybe $b is not on NCBI" if $?;

    if(!-e "$$settings{tempdir}/$b.fastq"){
      logmsg "Converting sra file to fastq with fastq-dump";
      system("fastq-dump -v -v -v --legacy-report -I --split-files -O $$settings{tempdir}/$b.fastq.tmp $$settings{tempdir}/$b.sra && mv $$settings{tempdir}/$b.fastq.tmp $$settings{tempdir}/$b.fastq");
      die "ERROR with fastq-dump" if $?;
    }

    if(!-e "$project/reads/$b.fastq.gz"){
      logmsg "Shuffling the reads";
      system("run_assembly_shuffleReads.pl $$settings{tempdir}/$b.fastq/*_1.fastq $$settings{tempdir}/$b.fastq/*_2.fastq | gzip -c > $project/reads/$b.fastq.gz.tmp && mv -v $project/reads/$b.fastq.gz.tmp $project/reads/$b.fastq.gz");
      die "ERROR with either run_assembly_shuffleReads.pl, gzip, or mv" if $?;
    }

    system("rm -rvf $$settings{tempdir}/$b.fastq/ $$settings{tempdir}/$b.sra");
    die "ERROR in removing these files" if $?;
  }
  return $newPath;
}
sub addAssembly{
  my($project,$settings)=@_;
  my $asm=File::Spec->rel2abs($$settings{'add-assembly'});
  if(-e $asm){
    my $symlink="$project/asm/".basename($asm);
    symlink($asm,$symlink);
    logmsg "$asm => $symlink";
  } else {
    logmsg "NOTE: I could not find file $asm";
    logmsg "I will see if it's on NCBI";

    my $b=basename($asm,qw(.fasta .fna .fa));
    my $fasta=_downloadAssembly($b,$project,$settings);
    move($fasta,"$project/asm/$b.fasta");
    die "Could not move file: $!" if $?;

    logmsg "Downloaded $b into $project/asm/$b.fasta";
  }
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
  if(-e $ref){
    copy($ref,$newPath) or die "ERROR: could not copy $ref to $newPath\n  $!";
    logmsg "cp $ref => $newPath";
  } else {
    my $b=basename($ref,qw(.fasta .fna .fa));
    my $fasta=_downloadAssembly($b,$project,$settings);
    move($fasta,$newPath);
    die "Could not move file: $!" if $?;

    logmsg "Downloaded $b into $newPath";
  }

  unlink("$project/reference/reference.fasta");
  symlink(basename($newPath),"$project/reference/reference.fasta");
  logmsg "symlink $newPath => $project/reference/reference.fasta";
  return 1;
}

sub _downloadAssembly{
  my($ID,$project,$settings)=@_;
  my $b=basename($ID,qw(.fasta .fna .fa));
  my $tmpFile="$project/$b.tmp.fasta";
  my $filteredFile="$project/$b.fasta";
  system("esearch -db nucleotide -query $b | efetch -format fasta > $tmpFile");
  die "ERROR with efetch or esearch" if $?;
  die "ERROR: downloaded a zero-byte file for $b. Is it the correct identifier?\n" if (-s "$project/$b.tmp.fasta" < 1);
  
  # filter the sequence
  my $in=Bio::SeqIO->new(-format=>"fasta",-file=>$tmpFile);
  my $out=Bio::SeqIO->new(-format=>"fasta",-file=>">$filteredFile");
  while(my $seq=$in->next_seq){
    next if($seq->length < 1);
    $out->write_seq($seq);
  }
  $in->close;
  $out->close;

  return $filteredFile;
}

sub deleteProject{
  my($project,$settings)=@_;
  for (@projectSubdir){
    system("rm -rfv '$project/$_' 1>&2"); die if $?;
  }
  system("rmdir -v $project 1>&2"); die if $?;
  return 1;
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

  "Manages a Lyve-SET project and confirms that a directory is a Lyve-SET project.
  Usage: $0 setProject [options]
  --create                        Create a set project
  --test DATASET                  Create a set project with test data (invokes --create). Options are:
                                    $testdata
  --delete                        Delete a set project
  --add-reads file.fastq.gz       Add reads to your project (using symlink)
  --add-reads SRR0123456          Add reads to your project (using NCBI/SRA/wget)
  --remove-reads file.fastq.gz    Remove reads from your project (using rm on reads, vcf, and bam directories)
  --add-assembly file.fasta       Add an assembly to your project (using symlink)
  --add-assembly NC_0123456       Add an assembly to your project (using NCBI/nucleotide/edirect)
  --remove-assembly file.fasta    Remove an assembly from your project (using rm)
  --change-reference file.fasta   Add a reference assembly to your project and remove any existing one (using cp)
  Examples
    set_manage.pl projDir --add-reads ../../reads/unknown/file.cleaned.fastq.gz
    set_manage.pl projDir --remove-reads file.cleaned.fastq.gz
  "
}
