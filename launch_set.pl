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
use Schedule::SGELK;

sub logmsg {local $0=basename $0;my $FH = *STDOUT; print $FH "$0: ".(caller(1))[3].": @_\n";}
my ($name,$scriptsdir,$suffix)=fileparse($0);
$scriptsdir=File::Spec->rel2abs($scriptsdir);

my $sge=Schedule::SGELK->new(-verbose=>1,-numnodes=>5,-numcpus=>8);
exit(main());

sub main{
  my $settings={trees=>1,clean=>1, msa=>1, mapper=>"smalt"};
  GetOptions($settings,qw(ref=s bamdir=s vcfdir=s tmpdir=s readsdir=s asmdir=s msadir=s help numcpus=s numnodes=i workingdir=s allowedFlanking=s keep min_alt_frac=s min_coverage=i trees! qsubxopts=s clean! msa! mapper=s snpcaller=s));
  $$settings{numcpus}||=8;
  $$settings{numnodes}||=6;
  $$settings{workingdir}||=$sge->get("workingdir");
  $$settings{allowedFlanking}||=0;
  $$settings{keep}||=0;
  $$settings{min_alt_frac}||=0.75;
  $$settings{min_coverage}||=10;
  $$settings{qsubxopts}||="";
  $$settings{mapper}=lc($$settings{mapper});
  $$settings{snpcaller}||="freebayes";

  logmsg "Checking to make sure all directories are in place";
  for my $param (qw(vcfdir bamdir msadir readsdir tmpdir asmdir)){
    my $b=$param;
    $b=~s/dir$//;
    $$settings{$param}||=$b;
    die "ERROR: Could not find $param under $$settings{$param}/ \n".usage() if(!-d $$settings{$param});
    $$settings{$param}=File::Spec->rel2abs($$settings{$param});
  }
  # SGE params
  for (qw(workingdir numnodes numcpus keep qsubxopts)){
    $sge->set($_,$$settings{$_});
  }

  die usage() if($$settings{help} || !defined($$settings{ref}) || !-f $$settings{ref});
  my $ref=$$settings{ref};

  logmsg "Simulating reads from any assemblies";
  simulateReads($settings);
  logmsg "Indexing the reference for read mapping";
  indexReference($ref,$settings);
  logmsg "Mapping reads";
  mapReads($ref,$$settings{readsdir},$$settings{bamdir},$settings);
  logmsg "Calling variants";
  variantCalls($ref,$$settings{bamdir},$$settings{vcfdir},$settings);

  if($$settings{msa}){
    logmsg "Creating a core hqSNP MSA";
    variantsToMSA($ref,$$settings{bamdir},$$settings{vcfdir},$$settings{msadir},$settings);
    logmsg "Launching set_processMsa.pl";
    $sge->pleaseExecute_andWait("set_processMsa.pl --auto --msaDir '$$settings{msadir}' --numcpus $$settings{numcpus}",{numcpus=>$$settings{numcpus},jobname=>"set_processMsa.pl"});
    #msaToPhylogeny($$settings{msadir},$settings) if($$settings{trees});
  }

  logmsg "Done!";

  return 0;
}

sub simulateReads{
  my($settings)=@_;
  my @asm=glob("$$settings{asmdir}/*.{fasta,fna,ffn,fa,mfa,mfasta}");
  return if(!@asm);

  my $exec=`which wombac-shred_fasta`;
  if($?){
    $ENV{PATH}="$ENV{PATH}:$FindBin::RealBin/lib/wombac-1.2/bin";
    $exec=`which wombac-shred_fasta`;
    die "ERROR could not find wombac-shred_fasta in your path" if $?;
  }
  chomp($exec);

  logmsg "Using Wombac to shred reads. Please give credit to Wombac if using this stage of SET. Details at http://www.vicbioinformatics.com/";
  for my $asm (@asm){
    my $b=fileparse $asm;
    my $outFastq="$$settings{readsdir}/$b.shredded.fastq";
    my $outGz="$outFastq.gz";
    next if(-e $outGz);
    logmsg "I did not find $outFastq simulated reads. Simulating now with $exec.";
    $sge->set("jobname","sim$b");
    $sge->pleaseExecute_andWait("$exec --verbose --readlen 1000 --coverage 50 --fastq --qual 40 '$asm' 2>&1 1> '$outFastq'");
    $sge->pleaseExecute_andWait("gzip -v '$outFastq' 2>&1");
  }
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

  if($$settings{mapper} eq 'smalt'){
    system("smalt index -k 5 -s 3 $ref $ref 2>&1");
    die if $?;
  } elsif($$settings{mapper} eq 'snap'){
    system("snap index $ref $ref.snap -s 16");
    die if $?;
  }
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
    my $clean=($$settings{clean})?"--clean":"--noclean"; # the clean parameter or not

    if($$settings{mapper} eq 'smalt'){
      $sge->pleaseExecute("$scriptsdir/launch_smalt.pl -ref $ref -f $fastq -b $bamPrefix.sorted.bam -tempdir $tmpdir --numcpus $$settings{numcpus} $clean",{jobname=>"smalt$b"});
    } elsif($$settings{mapper} eq 'snap'){
      $sge->pleaseExecute("$scriptsdir/launch_snap.pl -ref $ref -f $fastq -b $bamPrefix.sorted.bam -tempdir $tmpdir --numcpus $$settings{numcpus} $clean",{jobname=>"snap$b"});
    } else {
      die "ERROR: I do not understand the mapper $$settings{mapper}";
    }
  }
  logmsg "All mapping jobs have been submitted. Waiting on them to finish.";
  $sge->wrapItUp();
  return 1;
}

sub variantCalls{
  my($ref,$bamdir,$vcfdir,$settings)=@_;
  my @bam=glob("$bamdir/*.sorted.bam");
  my @jobid;
  for my $bam(@bam){
    my $b=fileparse($bam,".sorted.bam");
    $sge->set("jobname","varcall$b");
    if(-e "$vcfdir/$b.vcf"){
      logmsg "Found $vcfdir/$b.vcf. Skipping";
      next;
    }
    my $j; # job identifier
    if($$settings{snpcaller} eq 'freebayes'){
      $j=$sge->pleaseExecute("$scriptsdir/launch_freebayes.sh $ref $bam $vcfdir/$b.vcf $$settings{min_alt_frac} $$settings{min_coverage}",{numcpus=>1});
      # terminate called after throwing an instance of 'std::out_of_range'
    } elsif($$settings{snpcaller} eq 'callsam'){
      my $jobname="callsam$b";
      # call snps
      $j=$sge->pleaseExecute("$scriptsdir/lib/callsam/bin/callsam.pl $bam --numcpus $$settings{numcpus} --min-coverage $$settings{min_coverage} --min-frequency $$settings{min_alt_frac} --reference '$ref' > $vcfdir/unfiltered/$b.vcf",{numcpus=>$$settings{numcpus},jobname=>$jobname,qsubxopts=>""});
      # filter the vcf but make it depend on the callsam script to finish
      $sge->pleaseExecute("$scriptsdir/filterVcf.pl $vcfdir/unfiltered/$b.vcf --noindels -d 10 -o $vcfdir/$b.vcf",{qsubxopts=>"-hold_jid $jobname",numcpus=>1,jobname=>"filter$b"});
    }
    push(@jobid,$j);
  }
  logmsg "All variant-calling jobs have been submitted. Waiting on them to finish";
  $sge->wrapItUp();
  return 1;
}

sub variantsToMSA{
  my ($ref,$bamdir,$vcfdir,$msadir,$settings)=@_;
  my $logdir=$$settings{logdir};
  if(-e "$msadir/out.aln.fas.phy"){
    logmsg "Found $msadir/out.aln.fas.phy already present. Not re-converting.";
    return 1;
  }

  # Find the distance that we'd expect SNP-linkage, so that they can be filtered out
  if($$settings{allowedFlanking} eq 'auto'){
    my $allowedFlanking=`snpDistribution.pl $vcfdir/unfiltered/*.vcf`;
    die if $?;
    chomp($allowedFlanking);
    $$settings{allowedFlanking}=$allowedFlanking;
  }

  # find all "bad" sites
  my $bad="$vcfdir/allsites.txt";
  system("sort $vcfdir/*.badsites.txt | uniq > $bad"); die if $?;

  # convert VCFs to an MSA (long step)
  $sge->set("jobname","variantsToMSA");
  $sge->set("numcpus",$$settings{numcpus});
  $sge->pleaseExecute("vcfToAlignment.pl $bamdir/*.sorted.bam $vcfdir/*.vcf -o $msadir/out.aln.fas -r $ref -b $bad -a $$settings{allowedFlanking}");
  # convert VCFs to an MSA using a low-memory script
  $sge->pleaseExecute("vcfToAlignment_lowmem.pl $vcfdir/unfiltered/*.vcf $bamdir/*.sorted.bam -n $$settings{numcpus} -ref $ref -p $msadir/out_lowmem.aln.fas.pos.txt -t $msadir/out_lowmem.aln.fas.pos.tsv > $msadir/out_lowmem.aln.fas",{numcpus=>$$settings{numcpus},jobname=>"variantsToMSA_lowmem"});
  $sge->wrapItUp();

  # convert fasta to phylip and remove uninformative sites
  #$sge->set("jobname","msaToPhylip");
  #$sge->pleaseExecute_andWait("convertAlignment.pl -i $msadir/out.aln.fas -o $msadir/out.aln.fas.phy -f phylip -r");
  return 1;
}

sub msaToPhylogeny{
  my ($msadir,$settings)=@_;

  $sge->set("numcpus",$$settings{numcpus});
  # raxml: remove the previous run, but not if it finished successfully
  if(-e "$msadir/RAxML_info.out"){
    if(!-e "$msadir/RAxML_bipartitions.out"){
      unlink("$msadir/RAxML_info.out");
    }
  }
  if(!-e "$msadir/RAxML_info.out"){
    $sge->set("jobname","SET_raxml");
    my $rand =int(rand(999999999));
    my $rand2=int(rand(999999999));
    #$sge->pleaseExecute("(cd $msadir; raxmlHPC-PTHREADS -f a -s $msadir/out.aln.fas.phy -n out -T $$settings{numcpus} -m GTRGAMMA -N 100 -p $rand -x $rand2)");
    $sge->pleaseExecute("(cd $msadir; launch_raxml.sh out.aln.fas out)");
  }

  # phyml
  $sge->pleaseExecute("launch_phyml.sh $msadir/out.aln.fas.phy",{jobname=>"SET_phyml"});
  $sge->wrapItUp();
  return 1;
}

sub usage{
  $0=fileparse $0;
  "Usage: $0 -ref reference.fasta [-b bam/ -v vcf/ -t tmp/ -reads reads/ -m msa/ -asm asm/]
    Where parameters with a / are directories
    -r where fastq and fastq.gz files are located
    -b where to put bams
    -v where to put vcfs
    --msadir multiple sequence alignment and tree files (final output)
    -numcpus number of cpus
    -numnodes maximum number of nodes
    -w working directory where qsub commands can be stored. Default: CWD
    -all allowed flanking distance in bp. Nucleotides this close together cannot be considered as high-quality.
      Set -all to 'auto' to let SET determine this distance using snpDistribution.pl
    -asm directory of assemblies. Copy or symlink the reference genome assembly to use it if it is not already in the raw reads directory
    --nomsa to not make a multiple sequence alignment
    --notrees to not make phylogenies
    -q '-q long.q' extra options to pass to qsub. This is not sanitized.
    --noclean to not clean reads before mapping (faster, but you need to have clean reads to start with)
    --mapper smalt Choose the mapper you want to use. Choices: smalt (default), snap
    --snpcaller freebayes Choose the snp caller you want to use. Choices: freebayes (default), callsam
  "
}
