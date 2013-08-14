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
  GetOptions($settings,qw(ref=s bamdir=s vcfdir=s tmpdir=s readsdir=s msadir=s help numcpus=s numnodes=i workingdir=s allowedFlanking=i keep));
  $$settings{numcpus}||=8;
  $$settings{numnodes}||=6;
  $$settings{workingdir}||=$sge->get("workingdir");
  $$settings{allowedFlanking}||=0;
  $$settings{keep}||=0;

  logmsg "Checking to make sure all directories are in place";
  for my $param (qw(vcfdir bamdir msadir readsdir tmpdir)){
    my $b=$param;
    $b=~s/dir$//;
    $$settings{$param}||=$b;
    die "ERROR: Could not find $param under $$settings{$param}/ \n".usage() if(!-d $$settings{$param});
    $$settings{$param}=File::Spec->rel2abs($$settings{$param});
  }
  # SGE params
  for (qw(workingdir numnodes numcpus keep)){
    $sge->set($_,$$settings{$_});
  }

  die usage() if($$settings{help} || !defined($$settings{ref}) || !-f $$settings{ref});
  my $ref=$$settings{ref};

  indexReference($ref,$settings);
  logmsg "Mapping reads";
  mapReads($ref,$$settings{readsdir},$$settings{bamdir},$settings);
  logmsg "Calling variants";
  variantCalls($ref,$$settings{bamdir},$$settings{vcfdir},$settings);
  logmsg "Creating a core hqSNP MSA";
  variantsToMSA($ref,$$settings{bamdir},$$settings{vcfdir},$$settings{msadir},$settings);
  logmsg "MSA => phylogeny";
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
    my $bamPrefix="$bamdir/$b-".basename($ref,qw(.fasta .fna .fa));
    if(-e "$bamPrefix.sorted.bam"){
      logmsg "Found $bamPrefix.sorted.bam. Skipping.";
      next;
    }else{
      logmsg "Mapping to create $bamPrefix.sorted.bam";
    }
    $sge->set("jobname","map$b");
    $sge->pleaseExecute("$scriptsdir/launch_smalt.sh $ref $fastq $bamPrefix.sorted.bam $tmpdir");
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

  # find all "bad" sites
  my $bad="$vcfdir/allsites.txt";
  system("sort $vcfdir/*.badsites.txt | uniq > $bad"); die if $?;

  # convert VCFs to an MSA (long step)
  $sge->set("jobname","variantsToMSA");
  $sge->set("numcpus",$$settings{numcpus});
  $sge->pleaseExecute_andWait("vcfToAlignment.pl $bamdir/*.sorted.bam $vcfdir/*.vcf -o $msadir/out.aln.fas -r $ref -b $bad -a $$settings{allowedFlanking}");
  # convert fasta to phylip and remove uninformative sites
  $sge->set("jobname","msaToPhylip");
  $sge->pleaseExecute_andWait("convertAlignment.pl -i $msadir/out.aln.fas -o $msadir/out.aln.fas.phy -f phylip -r");
  return 1;
}

sub msaToPhylogeny{
  my ($msadir,$settings)=@_;
  $sge->set("numcpus",$$settings{numcpus});
  $sge->set("jobname","SET_raxml");
  my $rand =int(rand(999999999));
  my $rand2=int(rand(999999999));
  $sge->pleaseExecute("(cd $msadir; raxmlHPC-PTHREADS -f a -s $msadir/out.aln.fas.phy -n out -T $$settings{numcpus} -m GTRGAMMA -p $rand -x $rand2 -N 100)");
  $sge->set("jobname","SET_phyml");
  $sge->pleaseExecute("PhyML -i $msadir/out.aln.fas.phy -b -4 -m GTR -s BEST --quiet");
  $sge->wrapItUp();
  return 1;
}

sub usage{
  $0=fileparse $0;
  "Usage: $0 -ref reference.fasta [-b bam/ -v vcf/ -t tmp/ -reads reads/ -m msa/]
    Where parameters with a / are directories
    -r where fastq and fastq.gz files are located
    -b where to put bams
    -v where to put vcfs
    -m multiple sequence alignment and tree files (final output)
    -numcpus number of cpus
    -numnodes maximum number of nodes
    -w working directory where qsub commands can be stored. Default: CWD
    -a allowed flanking distance in bp. Nucleotides this close together cannot be considered as high-quality.
  "
}
