#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Spec;
use threads;
use Thread::Queue;

sub logmsg {local $0=basename $0;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}
my ($name,$scriptsdir,$suffix)=fileparse($0);
$scriptsdir=File::Spec->rel2abs($scriptsdir);

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(ref=s bamdir=s vcfdir=s tmpdir=s logdir=s readsdir=s msadir=s progressdir=s help numcpus=s numnodes=i));
  $$settings{numcpus}||=8;
  $$settings{numnodes}||=6;

  logmsg "Checking to make sure all directories are in place";
  for my $param (qw(logdir vcfdir bamdir msadir readsdir tmpdir progressdir)){
    my $b=$param;
    $b=~s/dir$//;
    $$settings{$param}||=$b;
    die "ERROR: Could not find $param under $$settings{$param}/ \n".usage() if(!-d $$settings{$param});
    $$settings{$param}=File::Spec->rel2abs($$settings{$param});
  }

  die usage() if($$settings{help} || !defined($$settings{ref}) || !-f $$settings{ref});
  my $ref=$$settings{ref};
  # TODO check if ref is readable and a fasta file

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
  my $tmpdir=$$settings{tmpdir};
  my $log=$$settings{logdir};
  my $progress=$$settings{progressdir};
  my @file=(glob("$readsdir/*.fastq"),glob("$readsdir/*.fastq.gz"));
  my @jobid;
  for my $fastq(@file){
    my $b=fileparse $fastq;
    unlink("$progress/smalt.$b.done");
    $$settings{jobname}="map$b";
    my $jobid=pleaseExecute("$scriptsdir/launch_smalt.sh $ref $fastq $bamdir/$b.bam $tmpdir",$settings);
    die if $?;
    push(@jobid,$jobid);
    waitOnJobs(\@jobid,$settings);
  }

  $$settings{mustfinish}=1;
  waitOnJobs(\@jobid,$settings);
  $$settings{mustfinish}=0;

  return 1;
}

sub variantCalls{
  my($ref,$bamdir,$vcfdir,$settings)=@_;
  my $progress=$$settings{progressdir};
  my @bam=glob("$bamdir/*.sorted.bam");
  my @jobid;
  my @donefile;
  #logmsg "DEBUG"; return 1;
  for my $bam(@bam){
    my $b=fileparse($bam,".sorted.bam");
    $$settings{jobname}="varcall$b";
    my $jobid=pleaseExecute("$scriptsdir/launch_freebayes.sh $ref $bam vcf/$b.vcf",$settings);
    push(@jobid,$jobid);
    #system("qsub -N 'q$b' -cwd -V -o log/$b.out -e log/$b.out $scriptsdir/launch_freebayes.sh $ref $bam vcf/$b.vcf");
    unlink("$progress/freebayes.$b.done");
    push(@donefile,"$progress/freebayes.$b.done");
    waitOnJobs(\@jobid,$settings);
  }
  $$settings{mustfinish}=1;
  waitOnJobs(\@jobid,$settings);
  $$settings{mustfinish}=0;
  
  # wait for variants to be finished
  die "There wasn't anything mapped!" if(!@donefile);
  for my $donefile (@donefile){
    logmsg "Waiting on $donefile";
    while(! -f $donefile){
      sleep 10;
      if(`qstat|wc -l` < 1){
        logmsg "WARNING: there is nothing left in the queue. Did the freebayes script die?";
      }
    }
  }
  return 1;
}

sub variantsToMSA{
  my ($ref,$bamdir,$vcfdir,$msadir,$settings)=@_;
  my $logdir=$$settings{logdir};
  my $progressdir=$$settings{progressdir};
  #logmsg "DEBUG"; return 1;
  system("qsub -pe smp $$settings{numcpus} -cwd -o $logdir/toMsa.log -e $logdir/toMsa.log -V $scriptsdir/launch_vcfToAlignment.sh $bamdir $vcfdir $ref $msadir/out.aln.fas");
  my $donefile="$progressdir/vcf2aln.done";
  unlink($donefile);
  logmsg "Variants are being transformed into an MSA. Waiting on $donefile to indicate that it is done.";
  while(! -f $donefile){
    sleep 2;
  }
  return 1;
}

sub msaToPhylogeny{
  my ($msadir,$settings)=@_;
  my $progressdir=$$settings{progressdir};
  my $donefile="$progressdir/raxml.done";
  system("cd $msadir; qsub -pe smp $$settings{numcpus} -cwd -o out -e out -V $scriptsdir/launch_raxml.sh out.aln.fas.phy out");
  # To get the job ID, parse the stdout of the qsub command which is something like
  #   Your job 123057 ("test.sh") has been submitted
  die if $?;
  while(! -f $donefile){
    sleep 2;
  }
  return 1;
}

# Pauses the program while we wait for these files.
# These files are usually zero-byte files to indicate a job finished.
sub waitOnDonefiles{
  my($donefile,$errfile,$settings)=@_;
  logmsg scalar(@$donefile)." in the queue." if(@$donefile >= 1 && @$donefile <=$$settings{numnodes});

  while(@$donefile > 1){
    for(my $i=0;$i<@$donefile;$i++){
      if(-e $$donefile[$i]){
        logmsg "A job finished: found $$donefile[$i]";
        splice(@$donefile,$i,1);
        last;
      }
    }
    sleep 1;
    last if(!$$settings{mustfinish} && @$donefile>$$settings{numnodes});
  }
  logmsg "DONE!" if($$settings{mustfinish});
  return @$donefile;
}

sub waitOnJobs{
  my($jobid,$settings)=@_;
  logmsg "We have reached node capacity ($$settings{numnodes})!" if(@$jobid >= $$settings{numnodes});
  while(@$jobid > 0){
    for(my $i=0;$i<@$jobid;$i++){
      my $state=checkJob($$jobid[$i],$settings);
      if($state eq 0){
        logmsg "A job finished: $$jobid[$i]";
        splice(@$jobid,$i,1);
        last;
      }
    }
    sleep 1;
    # break out if you don't have to finish yet but you can still add in another job
    last if(!$$settings{mustfinish} && @$jobid<$$settings{numnodes});
  }
  return @$jobid;
}

# submit to the cluster
  #system("qsub -pe smp $$settings{numcpus} -N 'q$b' -cwd -V -o $log/$b.out -e $log/$b.out $scriptsdir/launch_smalt.sh $ref $fastq $bamdir/$b.bam $tmpdir");
sub pleaseExecute{
  my($cmd,$settings)=@_;
  local $0=basename $0;
  my $jobid=-1; # -1 is an error state
  $$settings{jobname}||="run$0";
  $$settings{logfile}||="$0.log";
  $$settings{numcpus}||=1;
  $cmd="echo '$cmd' | qsub -pe smp $$settings{numcpus} -N $$settings{jobname} -cwd -V -o $$settings{logfile} -e $$settings{logfile} 2>&1";
  my $out=`$cmd`;
  print $out;
  if($out=~/Your job (\d+)/){
    $jobid=$1;
  } else {
    logmsg "WARNING: the last job submitted did not have an obvious jobid. It can't be tracked!";
  }
  return $jobid;
}

# return the job status of the id
sub checkJob{
  my($jobid,$settings)=@_;
  my $state=0;
  open(QSTAT,"qstat|") or die "ERROR: could not execute qstat!";
  while(<QSTAT>){
    s/^\s+|\s+$//g;
    my @F=split /\s+/;
    if($F[0] eq $jobid){
      $state=$F[4];
    }
  }
  close QSTAT;
  return $state;
}

sub usage{
  $0=fileparse $0;
  "Usage: $0 -ref reference.fasta [-b bam/ -v vcf/ -t tmp/ -l log/ -reads reads/ -m msa/ -p progress/]
    Where parameters with a / are directories
    -r where fastq and fastq.gz files are located
    -b where to put bams
    -v where to put vcfs
    -m multiple sequence alignment and tree files (final output)
    -numcpus number of cpus (default: 8)
    -numnodes maximum number of nodes (default: 6)
  "
}
