#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::RealBin/../lib";

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
#use ExtUtils::Command;


my ($name,$scriptsdir,$suffix)=fileparse($0);
$scriptsdir=File::Spec->rel2abs($scriptsdir);

# Logging
my $logmsgFh;
sub logmsg {
  local $0=basename $0;
  my $FH = *STDOUT; 
  my $msg="$0: ".(caller(1))[3].": @_\n";

  # print the message to the logfile if it's not the same as stdout
  #print join("\t","===",fileno($logmsgFh),fileno($FH))."\n";
  if(defined($logmsgFh) && fileno($logmsgFh) > 0){
    print $logmsgFh $msg;
  }
  print $FH $msg;
}
local $SIG{'__DIE__'} = sub { local $0=basename $0; my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; my $msg="$0: ".(caller(1))[3].": ".$e; logmsg($msg);  die(); };
END{
  close $logmsgFh if(defined($logmsgFh) && fileno($logmsgFh) > 0);
}

my $sge=Schedule::SGELK->new(-verbose=>1,-numnodes=>20,-numcpus=>8);
exit(main());

sub main{
  # start with the settings that are on by default, and which can be turned off by, e.g., --noclean
  my $settings={trees=>1,msa=>1, matrix=>1, clean=>1};
  GetOptions($settings,qw(ref=s bamdir=s logdir=s vcfdir=s tmpdir=s readsdir=s asmdir=s msadir=s help numcpus=s numnodes=i allowedFlanking=s keep min_alt_frac=s min_coverage=i trees! queue=s qsubxopts=s msa! matrix! mapper=s snpcaller=s msa-creation=s clean)) or die $!;
  # Lyve-SET
  $$settings{allowedFlanking}||=0;
  $$settings{keep}||=0;
  $$settings{min_alt_frac}||=0.75;
  $$settings{min_coverage}||=10;
  # modules defaults
  $$settings{mapper}||="smalt";
  $$settings{snpcaller}||="freebayes";
  $$settings{'msa-creation'}||="lyve-set-lowmem";
  # queue stuff
  $$settings{vcfToAlignment_xopts}||="-l mem_free=100G -q highmem.q";
  $$settings{numcpus}||=1;
  $$settings{numnodes}||=20;
  $$settings{qsubxopts}||="";
  $$settings{queue}||="";
  # Some things need to just be lowercase to make things easier downstream
  $$settings{$_}=lc($$settings{$_}) for(qw(msa-creation snpcaller mapper));

  # A warning about --noclean
  if(!$$settings{clean}){
    die "Warning: the --noclean option has been removed in Lyve-SET v0.9.2.  You should clean your reads before using Lyve-SET.";
  }

  ##########################################################
  ### Other defaults: reference genome; default directories#
  ##########################################################
  my $project=shift(@ARGV) || '.'; # by default the user is in the project directory
  die usage($settings) if($$settings{help});
  logmsg "Project was defined as $project"; # No period at the end of sentence to avoid ambiguous '..' in the logmsg
  # Checks to make sure that this is a project
  #system("set_manage.pl '$project'"); die if $?;
  # Set the defaults if they aren't set already
  $$settings{ref}||="$project/reference/reference.fasta";
  # Check SET directories' existence and set their defaults
  for my $param (qw(vcfdir bamdir msadir readsdir tmpdir asmdir logdir)){
    my $b=$param;
    $b=~s/dir$//;  # e.g., vcfdir => vcf
    $$settings{$param}||="$project/$b";
    die "ERROR: Could not find $param under $$settings{$param}/\n  mkdir $$settings{$param} to resolve this problem.\n".usage($settings) if(!-d $$settings{$param});
    $$settings{$param}=File::Spec->rel2abs($$settings{$param});
  }

  # SGE params
  $$settings{workingdir}=$$settings{logdir};
  for (qw(workingdir numnodes numcpus keep qsubxopts queue)){
    $sge->set($_,$$settings{$_}) if($$settings{$_});
  }

  # open the log file now that we know the logdir is there
  open($logmsgFh,">$$settings{logdir}/launch_set.log") or die "ERROR: could not open log file $$settings{logdir}/launch_set.log";
  logmsg "======\nOpened logfile at $$settings{logdir}/launch_set.log\n=====";

  # Check the reference parameter
  die "ERROR: reference file was not given\n".usage($settings) if(!defined($$settings{ref}));
  die "ERROR: Could not find the reference file at $$settings{ref}\n".usage($settings) if(!-f $$settings{ref});
  my $ref=$$settings{ref};

  simulateReads($settings);
  indexReference($ref,$settings);
  mapReads($ref,$$settings{readsdir},$$settings{bamdir},$settings);
  variantCalls($ref,$$settings{bamdir},$$settings{vcfdir},$settings);

  if($$settings{matrix}){
    variantsToMatrix($ref,$$settings{bamdir},$$settings{vcfdir},$$settings{msadir},$settings);
  } else {
    logmsg "The matrix was not requested; wrapping up";
    return 0;
  }

  if($$settings{msa}){
    matrixToAlignment($settings);
  } else {
    logmsg "The alignment was not requested; wrapping up";
    return 0;
  }

  if($$settings{trees}){
    logmsg "Launching set_processMsa.pl";
    $sge->pleaseExecute("set_processMsa.pl --auto --msaDir '$$settings{msadir}' --numcpus $$settings{numcpus} 2>&1 | tee $$settings{logdir}/set_processMsa.log ",{numcpus=>$$settings{numcpus},jobname=>"set_processMsa.pl"});
    $sge->wrapItUp();
  } else {
    logmsg "The phylogeny was not requested; wrapping up";
    return 0;
  }

  return 0;
}

sub simulateReads{
  my($settings)=@_;
  logmsg "Simulating reads from any assemblies";
  my @asm=glob("$$settings{asmdir}/*.{fasta,fna,ffn,fa,mfa,mfasta}");
  return if(!@asm);

  my $exec=`which wgsim`;
  die "ERROR: Could not find wgsim which is part of samtools" if $?;
  chomp($exec);

  for my $asm (@asm){
    my $b=fileparse $asm;
    my $outFastq="$$settings{readsdir}/$b.wgsim.fastq";
    my $outGz="$outFastq.gz";
    my $out1="$$settings{tmpdir}/$b.1.fastq";
    my $out2="$$settings{tmpdir}/$b.2.fastq";
    next if(-e $outGz && !-e $outFastq); 
    logmsg "I did not find $outFastq simulated reads. Simulating now with $exec.";
    
    # Four successive jobs that depend on each other:
    #   1. Simulate into paired ends
    #   2. Shuffle
    #   3. Gzip
    #   4. Remove simulated unshuffled reads (depends on the shuffle step and not gzip)
    $sge->pleaseExecute("$exec -d 400 -N 100000 -1 250 -2 250 -e 0.02 -r 0.0009 -R 0.0001 -h '$asm' $out1 $out2",{jobname=>"wgsim$b",numcpus=>1});
    $sge->pleaseExecute("run_assembly_shuffleReads.pl $out1 $out2 > $outFastq",{jobname=>"shuffle$b",qsubxopts=>"-hold_jid wgsim$b",numcpus=>1});
    $sge->pleaseExecute("gzip -v '$outFastq' 2>&1",{jobname=>"gzip$b",qsubxopts=>"-hold_jid shuffle$b",numcpus=>1});
    $sge->pleaseExecute("rm -v $out1 $out2",{jobname=>"rmTmp$b",qsubxopts=>"-hold_jid shuffle$b",numcpus=>1});
  }
  $sge->wrapItUp();
}

sub indexReference{
  my($ref,$settings)=@_;
  logmsg "Indexing the reference for read mapping";

  # sanity check: see if the reference has dashes in its defline
  my $in=Bio::SeqIO->new(-file=>$ref);
  while(my $seq=$in->next_seq){
    my $defline=$seq->id." ".$seq->desc;
    die "Dashes are not allowed in the defline\n Offending defline: $defline" if($defline=~/\-/);
  }

  logmsg "Indexing with $$settings{mapper}";
  if($$settings{mapper} eq 'smalt'){
    return $ref if(-e "$ref.sma" && -e "$ref.smi");
    system("smalt index -k 5 -s 3 $ref $ref 2>&1");
    die if $?;
  } elsif($$settings{mapper} eq 'snap'){
    return $ref if(-d "$ref.snap" && -e "$ref.snap/GenomeIndex");
    system("snap index $ref $ref.snap -s 16");
    die if $?;
  }
  return $ref;
}

sub mapReads{
  my($ref,$readsdir,$bamdir,$settings)=@_;
  logmsg "Mapping reads with $$settings{mapper}";
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
      $sge->pleaseExecute("$scriptsdir/launch_smalt.pl -ref $ref -f $fastq -b $bamPrefix.sorted.bam -tempdir $tmpdir --numcpus $$settings{numcpus} ",{jobname=>"smalt$b"});
    } elsif($$settings{mapper} eq 'snap'){
      $sge->pleaseExecute("$scriptsdir/launch_snap.pl -ref $ref -f $fastq -b $bamPrefix.sorted.bam -tempdir $tmpdir --numcpus $$settings{numcpus} ",{jobname=>"snap$b"});
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
  logmsg "Calling variants with $$settings{snpcaller}";
  my @bam=glob("$bamdir/*.sorted.bam");

  # see if bgzip and tabix exist
  my $bgzip=`which bgzip 2>/dev/null`; chomp($bgzip);
  my $tabix=`which tabix 2>/dev/null`; chomp($tabix);

  # TODO make the variant callers output to bgzip and indexed files
  for my $bam(@bam){
    my $b=fileparse($bam,".sorted.bam");
    if(-e "$vcfdir/$b.vcf" || -e "$vcfdir/$b.vcf.gz"){
      logmsg "Found $vcfdir/$b.vcf. Skipping";
      next;
    }
    logmsg "Calling SNPs into $vcfdir/$b.vcf";
    my $jobname=""; # the job that the filter/sort scripts are waiting for
    if($$settings{snpcaller} eq 'freebayes'){
      $jobname="freebayes$b";
      $sge->pleaseExecute("$scriptsdir/launch_freebayes.sh $ref $bam $vcfdir/$b.vcf $$settings{min_alt_frac} $$settings{min_coverage}",{numcpus=>1,jobname=>$jobname});
      # terminate called after throwing an instance of 'std::out_of_range'
    } elsif($$settings{snpcaller} eq 'varscan'){
      $jobname="varscan$b";
      $sge->pleaseExecute("$scriptsdir/launch_varscan.pl $bam --reference $ref > $vcfdir/unfiltered/$b.vcf",{numcpus=>1,jobname=>$jobname,qsubxopts=>""});
      # sort VCF
      $sge->pleaseExecute("mv $vcfdir/unfiltered/$b.vcf $vcfdir/unfiltered/$b.vcf.tmp && vcf-sort < $vcfdir/unfiltered/$b.vcf.tmp > $vcfdir/unfiltered/$b.vcf",{jobname=>"sort$b",qsubxopts=>"-hold_jid $jobname",numcpus=>1});
      # filter VCF
      $sge->pleaseExecute("$scriptsdir/filterVcf.pl $vcfdir/unfiltered/$b.vcf --noindels -d $$settings{min_coverage} -o $vcfdir/$b.vcf",{qsubxopts=>"-hold_jid sort$b",numcpus=>1,jobname=>"filter$b"});
      $jobname="filter$b";
    } elsif($$settings{snpcaller} eq 'callsam'){
      $jobname="callsam$b";
      # call snps
      $sge->pleaseExecute("$scriptsdir/../lib/callsam/bin/callsam_MT.pl $bam --numcpus $$settings{numcpus} --min-coverage $$settings{min_coverage} --min-frequency $$settings{min_alt_frac} --reference '$ref' > $vcfdir/unfiltered/$b.vcf",{numcpus=>$$settings{numcpus},jobname=>$jobname,qsubxopts=>""});
      # sort VCF
      $sge->pleaseExecute("mv $vcfdir/unfiltered/$b.vcf $vcfdir/unfiltered/$b.vcf.tmp && vcf-sort < $vcfdir/unfiltered/$b.vcf.tmp > $vcfdir/unfiltered/$b.vcf",{jobname=>"sort$b",qsubxopts=>"-hold_jid $jobname",numcpus=>1});
      # filter VCF
      $jobname="filter$b";
      $sge->pleaseExecute("$scriptsdir/filterVcf.pl $vcfdir/unfiltered/$b.vcf --noindels -d $$settings{min_coverage} -o $vcfdir/$b.vcf",{qsubxopts=>"-hold_jid sort$b",numcpus=>1,jobname=>$jobname});
    } else {
      die "ERROR: I do not understand snpcaller $$settings{snpcaller}";
    }

    # TODO move filtering here

    # bgzip and tabix indexing
    # TODO enable this if statement or put it into the launch_* scripts
    if($bgzip && $tabix){
      indexAndCompressVcf("$vcfdir/unfiltered/$b.vcf",$jobname,$settings);
      indexAndCompressVcf("$vcfdir/$b.vcf",$jobname,$settings);
    }
    else {
      logmsg "Either bgzip or tabix is not in your path, and so I will not compress and index VCFs";
    }
  } # END bam while loop
  logmsg "All variant-calling jobs have been submitted. Waiting on them to finish";
  $sge->wrapItUp();
  return 1;
}

# returns the last SGE job information in a hash ref
sub indexAndCompressVcf{
  my($vcf,$holdjid,$settings)=@_;
  my $j={};
  eval{
    $j=$sge->pleaseExecute("
      vcf-sort < '$vcf' > '$vcf.sorted.tmp' && mv '$vcf.sorted.tmp' '$vcf' && \
      bgzip -f '$vcf' && \
      tabix '$vcf.gz'
    ",{qsubxopts=>"-hold_jid $holdjid",jobname=>"sortAndCompress",numcpus=>1});
  };
  if($@){
    logmsg "Warning: Could not compress and index $vcf: $@";
  }
  return $j;
}

sub variantsToMatrix{
  my ($ref,$bamdir,$vcfdir,$msadir,$settings)=@_;
  logmsg "Creating a core hqSNP matrix";
  my $logdir=$$settings{logdir};

  my $matrix="$msadir/out.matrix.tsv";
  my $informative="$msadir/informative.matrix.tsv";
  if(-e $informative){
    logmsg "Found $informative -- already present. Not re-converting.";
    return 1;
  }

  # input files
  my @vcf=glob("$vcfdir/*.vcf $vcfdir/unfiltered/*.vcf.gz");
  my @bam=glob("$bamdir/*.sorted.bam");
  my @infile=(@bam,@vcf);
  my $infile="'".join("' '",@infile)."'";

  # Find the distance that we'd expect SNP-linkage, so that they can be filtered out
  if($$settings{allowedFlanking} eq 'auto'){
    my $allowedFlanking=`snpDistribution.pl @vcf`;
    die if $?;
    chomp($allowedFlanking);

    # let's make it a little bit more strict actually.
    $$settings{allowedFlanking}=$allowedFlanking*3;
  }
  logmsg "AllowedFlanking: $$settings{allowedFlanking}";

  my $j;
  $j=$sge->pleaseExecute("vcfToMatrix.pl $infile -r $ref --numcpus $$settings{numcpus} > $matrix",{jobname=>"vcfToMatrix",numcpus=>$$settings{numcpus}});
  # shorten the deflines to remove the directory names
  $j=$sge->pleaseExecute("sed -i.bak 's|^.*/||g' '$matrix'",{qsubxopts=>"-hold_jid $$j{jobid}",jobname=>"sedOnAln",numcpus=>1});
  $j=$sge->pleaseExecute("removeUninformativeSitesFromMatrix.pl --min-distance $$settings{allowedFlanking} < $matrix > $informative",{qsubxopts=>"-hold_jid $$j{jobid}",numcpus=>1,jobname=>"removeUninformativeSitesFromMatrix"});

  $sge->wrapItUp();
  return $informative;
}

sub matrixToAlignment{
  my($settings)=@_;
  my $outMsa="$$settings{msadir}/out.aln.fas";
  my $matrix="$$settings{msadir}/informative.matrix.tsv";
  if(-e $outMsa){
    logmsg "Found $outMsa and so I will not remake it";
    return $outMsa;
  }

  $sge->pleaseExecute("matrixToAlignment.pl < $matrix > $outMsa",{jobname=>"matrixToAlignment",numcpus=>1});
  $sge->wrapItUp();

  return $outMsa;
}

sub usage{
  my($settings)=@_;

  ## Format a few variables correctly
  # simplified pathnames for some options
  my @dir=qw(asmdir msadir readsdir bamdir vcfdir tmpdir logdir);
  $$settings{$_}=File::Spec->abs2rel($_).'/' for(@dir);
  # right padding for some options
  $$settings{$_}=reverse(sprintf("%15s","".reverse($$settings{$_}))) for(qw(mapper snpcaller msa-creation allowedFlanking),@dir);
  #$$settings{$_}=sprintf("%10s",$$settings{$_}) for(qw(mapper snpcaller msa-creation));
  $0=fileparse $0;

  # The help menu
  my $help="$0: Launches the Lyve-SET pipeline
    Usage: $0 [project] [-ref reference.fasta]
    If project is not given, then it is assumed to be the current working directory.
    If reference is not given, then it is assumed to be proj/reference/reference.fasta
    Where parameters with a / are directories
    -ref      proj/reference/reference.fasta   The reference genome assembly
    -reads    $$settings{readsdir} where fastq and fastq.gz files are located
    -bam      $$settings{bamdir} where to put bams
    -vcf      $$settings{vcfdir} where to put vcfs
    --tmpdir  $$settings{tmpdir} tmp/ Where to put temporary files
    --msadir  $$settings{msadir} multiple sequence alignment and tree files (final output)
    --logdir  $$settings{logdir} Where to put log files. Qsub commands are also stored here.
    -asm      $$settings{asmdir} directory of assemblies. Copy or symlink the reference genome assembly to use it if it is not already in the raw reads directory

    SNP MATRIX OPTIONS
    --allowedFlanking  $$settings{allowedFlanking} allowed flanking distance in bp. Nucleotides this close together cannot be considered as high-quality.  Set to -1 to let SET determine this distance using snpDistribution.pl
    --min_alt_frac     $$settings{min_alt_frac}  The percent consensus that needs to be reached before a SNP is called. Otherwise, 'N'
    --min_coverage     $$settings{min_coverage}  Minimum coverage needed before a SNP is called. Otherwise, 'N'
    ";
    return "$help\n  --help To view more help\n" if(!$$settings{help});

    $help.="
    SKIP CERTAIN STEPS
    --nomatrix to not create an hqSNP matrix
    --nomsa to not make a multiple sequence alignment
    --notrees to not make phylogenies
    MODULES
    --mapper       $$settings{mapper}   Which mapper? Choices: smalt, snap
    --snpcaller    $$settings{snpcaller}   Which SNP caller? Choices: freebayes, callsam, varscan
    SCHEDULER AND MULTITHREADING OPTIONS
    --queue     all.q         The default queue to use.
    --qsubxopts '-N lyve-set' extra options to pass to qsub. This is not sanitized; internal options might overwrite yours.
    --numnodes  20  maximum number of nodes
    --numcpus   $$settings{numcpus}  number of cpus
  ";
  return $help;
}
