#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::RealBin/lib";

use strict;
use warnings;
use Bio::Perl;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Spec;
use threads;
use Thread::Queue;
#use Schedule::SGE;
use Schedule::SGELK;

sub logmsg {local $0=basename $0;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}
my ($name,$scriptsdir,$suffix)=fileparse($0);
$scriptsdir=File::Spec->rel2abs($scriptsdir);

my $sge=Schedule::SGELK->new(-verbose=>1,-numnodes=>5,-numcpus=>8);
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(ref=s bamdir=s vcfdir=s tmpdir=s logdir=s readsdir=s msadir=s help numcpus=s numnodes=i workingdir=s));
  $$settings{numcpus}||=8;
  $$settings{numnodes}||=6;
  $$settings{workingdir}||=$sge->get("workingdir");

  logmsg "Checking to make sure all directories are in place";
  for my $param (qw(logdir vcfdir bamdir msadir readsdir tmpdir)){
    my $b=$param;
    $b=~s/dir$//;
    $$settings{$param}||=$b;
    die "ERROR: Could not find $param under $$settings{$param}/ \n".usage() if(!-d $$settings{$param});
    $$settings{$param}=File::Spec->rel2abs($$settings{$param});
  }
  # SGE params
  for (qw(workingdir numnodes numcpus)){
    $sge->set($_,$$settings{$_});
  }

  die usage() if($$settings{help} || !defined($$settings{ref}) || !-f $$settings{ref});
  my $ref=$$settings{ref};

  indexReference($ref,$settings);
  logmsg "Mapping reads";
  mapReads($ref,$$settings{readsdir},$$settings{bamdir},$settings);
  variantCalls($ref,$$settings{bamdir},$$settings{vcfdir},$settings);
  variantsToMSA($ref,$$settings{bamdir},$$settings{vcfdir},$$settings{msadir},$settings);
  msaToPhylogeny($$settings{msadir},$settings);

  logmsg "Done!";

  return 0;
}

sub indexReference{
  my($ref,$settings)=@_;

  return $ref if(-e "$ref.sma" && -e "$ref.smi");
  # sanity check: see if the reference has dashes in its defline
  my $in=Bio::SeqIO->new(-file=>$ref);
  while(my $seq=$in->next_seq){
    my $defline=$seq->id." ".$seq->desc;
    die "Dashes are not allowed in the defline\n Offending defline: $defline" if($defline=~/\-/);
  }
  system("smalt index -k 5 -s 3 $ref $ref 2>&1");
  die if $?;
  return $ref;
}

sub mapReads{
  my($ref,$readsdir,$bamdir,$settings)=@_;
  $sge->set("numcpus",$$settings{numcpus});
  my $tmpdir=$$settings{tmpdir};
  my $log=$$settings{logdir};
  my @file=(glob("$readsdir/*.fastq"),glob("$readsdir/*.fastq.gz"));
  my @job;
  for my $fastq(@file){
    my $b=fileparse $fastq;
    $sge->set("jobname","map$b");
    $sge->pleaseExecute("$scriptsdir/launch_smalt.sh $ref $fastq $bamdir/$b.bam $tmpdir");
  }
  logmsg "All mapping jobs have been submitted. Waiting on them to finish.";
  $sge->wrapItUp();
  return 1;
}

sub variantCalls{
  my($ref,$bamdir,$vcfdir,$settings)=@_;
  $sge->set("numcpus",1);
  my @bam=glob("$bamdir/*.sorted.bam");
  my @jobid;
  for my $bam(@bam){
    my $b=fileparse($bam,".sorted.bam");
    $sge->set("jobname","varcall$b");
    $sge->pleaseExecute("$scriptsdir/launch_freebayes.sh $ref $bam vcf/$b.vcf");
    #system("qsub -N 'q$b' -cwd -V -o log/$b.out -e log/$b.out $scriptsdir/launch_freebayes.sh $ref $bam vcf/$b.vcf");
  }
  logmsg "All variant-calling jobs have been submitted. Waiting on them to finish";
  $sge->wrapItUp();
  return 1;
}

sub variantsToMSA{
  my ($ref,$bamdir,$vcfdir,$msadir,$settings)=@_;
  my $logdir=$$settings{logdir};
  $sge->set("jobname","variantsToMSA");
  $sge->set("numcpus",$$settings{numcpus});
  logmsg "Creating a core hqSNP MSA";
  $sge->pleaseExecute_andWait("$scriptsdir/launch_vcfToAlignment.sh $bamdir $vcfdir $ref $msadir/out.aln.fas");
  #system("qsub -pe smp $$settings{numcpus} -cwd -o $logdir/toMsa.log -e $logdir/toMsa.log -V $scriptsdir/launch_vcfToAlignment.sh $bamdir $vcfdir $ref $msadir/out.aln.fas");
  return 1;
}

sub msaToPhylogeny{
  my ($msadir,$settings)=@_;
  $sge->set("numcpus",$$settings{numcpus});
  $sge->set("jobname","msaToPhylogeny");
  logmsg "Inferring the phylogeny from the MSA";
  $sge->pleaseExecute_andWait("(cd $msadir; $scriptsdir/launch_raxml.sh out.aln.fas.phy out)");
  #system("cd $msadir; qsub -pe smp $$settings{numcpus} -cwd -o out -e out -V $scriptsdir/launch_raxml.sh out.aln.fas.phy out");
  return 1;
}

sub usage{
  $0=fileparse $0;
  "Usage: $0 -ref reference.fasta [-b bam/ -v vcf/ -t tmp/ -l log/ -reads reads/ -m msa/]
    Where parameters with a / are directories
    -r where fastq and fastq.gz files are located
    -b where to put bams
    -v where to put vcfs
    -m multiple sequence alignment and tree files (final output)
    -numcpus number of cpus
    -numnodes maximum number of nodes
    -w working directory where qsub commands can be stored. Default: CWD
  "
}
