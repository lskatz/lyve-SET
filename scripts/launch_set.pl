#!/usr/bin/env perl

require 5.12.0;  # includes the "not implemented/yadda yadda operator (...)"

# Update paths
use FindBin qw/$RealBin/;
use lib "$FindBin::RealBin/../lib";
# Make sure the path is set up correctly but do not take priority 
# over what the user truly wants, as dictated by the already-existing path.
# This cannot be used as a solution if the user has the incorrect
# software version in the path (e.g., samtools 0.1.19)
#$ENV{PATH}="$ENV{PATH}:$FindBin::RealBin";
# Aw, screw it.  Force the path on the user.
$ENV{PATH}="$RealBin:$ENV{PATH}";

use lib "$FindBin::RealBin/../lib/lib/perl5";

use strict;
use warnings;
use POSIX qw/strftime/;
use Bio::Perl;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw/rel2abs abs2rel/;
use File::Temp qw/tempdir/;
use threads;
use Thread::Queue;
use Schedule::SGELK;
use Config::Simple;
use Array::IntSpan;

use LyveSET qw/@fastaExt @richseqExt @fastqExt @bamExt @vcfExt rangeUnion rangeInversion/;

my ($name,$scriptsdir,$suffix)=fileparse($0);
$scriptsdir=rel2abs($scriptsdir);

# Logging
my $logmsgFh;
sub logmsg {
  local $0=basename $0;
  my $FH = *STDOUT; 
  my $caller=(caller(1))[3];
  $caller=(caller(2))[3] if($caller=~/__ANON__/); # Don't care about __ANON__ subroutines
  $caller=~s/(main:*)+//;
  #$caller=~s/^__ANON__//;  # Don't care about __ANON__ subroutines
  $caller=~s/^\s+|\s+$//g; # trim
  $caller="$caller:" if($caller ne "");

  my $msg="$0: $caller @_\n";

  # print the message to the logfile if it's not the same as stdout
  #print join("\t","===",fileno($logmsgFh),fileno($FH))."\n";
  if(defined($logmsgFh) && fileno($logmsgFh) > 0){
    print $logmsgFh $msg;
  }
  print $FH $msg;
}
#local $SIG{'__DIE__'} = sub { local $0=basename $0; my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; logmsg($e); die(); };
END{
  close $logmsgFh if(defined($logmsgFh) && fileno($logmsgFh) > 0);
}

my $invocation="$0 ".join(" ",@ARGV);
my $sge=Schedule::SGELK->new(-verbose=>1,-numnodes=>20,-numcpus=>1);
exit(main());

sub main{
  # Initialize settings by reading the global configuration
  # file, LyveSET.conf
  my $settings=readGlobalSettings(1);
  GetOptions($settings,qw(ref=s bamdir=s logdir=s vcfdir=s tmpdir=s readsdir=s asmdir=s msadir=s help numcpus=s numnodes=i allowedFlanking=s keep min_alt_frac=s min_coverage=i trees! queue=s qsubxopts=s msa! matrix! mapper=s snpcaller=s mask-phages! mask-cliffs! fast downsample sample-sites singleend presets=s read_cleaner=s qsub!)) or die $!;

  # What options change when --fast is used?
  if($$settings{fast}){
    $$settings{mapper}="snap";
    $$settings{'mask-phages'}=0;
    $$settings{'mask-cliffs'}=0;
    $$settings{downsample}=1;
    $$settings{'sample-sites'}=1;
  }

  # If the user wants to use preset configurations
  if($$settings{presets}){
    # Load the configuration
    my $confPath="$FindBin::RealBin/../config/presets.conf";
    my $cfg=new Config::Simple($confPath) or die "ERROR: could not open $confPath";
    # Find the configuration that was specified
    my $block=$cfg->block($$settings{presets});
    die "ERROR: configuration $$settings{presets} was not found in $confPath!" if(!keys(%$block));
    # Put the presets into the current settings
    $$settings{$_}=$$block{$_} for(keys(%$block));
  }

  # Some things need to just be lowercase to make things easier downstream
  $$settings{$_}=lc($$settings{$_}) for(qw(snpcaller mapper read_cleaner));

  ###########################################################
  ### Other defaults: reference genome; default directories #
  ###########################################################
  my $project=shift(@ARGV) || '.'; # by default the user is in the project directory
  die usage($settings) if($$settings{help});
  logmsg "Project was defined as $project";
  # Checks to make sure that this is a project
  #system("set_manage.pl '$project'"); die if $?;

  # Set the defaults if they aren't set already
  #$$settings{ref}||="$project/reference/reference.fasta";
  if($$settings{ref}){
    $$settings{ref}=$$settings{ref}; # if set, then leave it alone.
  } elsif(-e "$project/reference/reference.gbk"){
    $$settings{ref}="$project/reference/reference.gbk.fasta"; # this file might not exist yet and will be created soon
  } elsif(-e "$project/reference/reference.fasta"){
    $$settings{ref}="$project/reference/reference.fasta";
  } else {
    die "ERROR: reference file was not given\n".usage($settings);
  }
  #
  # If a richseq (gbk/embl) were given as the reference, convert it to fasta and update the variable to the fasta
  $$settings{richseq}||=""; # the genbank or embl file
  my($refname,$refdir,$refext)=fileparse($$settings{ref},@richseqExt);
  if($refext=~/gbk$|gb$|embl$/){
    $$settings{richseq}=$$settings{ref};
    $$settings{ref}=$$settings{richseq}.".fasta";
    die "ERROR: $$settings{richseq} does not exist" if(!-e $$settings{richseq});
    richseqToFasta($$settings{richseq},$$settings{ref},$settings);
  }

  # Ok, it's a defined string, but does the file exist?
  die "ERROR: Could not find the reference file at $$settings{ref}\n".usage($settings) if(!-e $$settings{ref});
  my $ref=$$settings{ref}; # $ref is called so many times, it might as well be its own variable
  $$settings{refdir}||=dirname($$settings{ref});

  # Check SET directories' existence and set their defaults
  die "ERROR: could not find dir $project" if(!-e $project);
  for my $param (qw(vcfdir bamdir msadir readsdir tmpdir asmdir logdir)){
    my $b=$param;
    $b=~s/dir$//;  # e.g., vcfdir => vcf
    $$settings{$param}||="$project/$b";
    die "ERROR: Could not find $param under $$settings{$param}/\n  mkdir $$settings{$param} to resolve this problem.\n".usage($settings) if(!-d $$settings{$param});
    $$settings{$param}=rel2abs($$settings{$param});
  }

  # SGE params
  mkdir "$$settings{logdir}/SGELK" if(!-d "$$settings{logdir}/SGELK");
  $$settings{workingdir}="$$settings{logdir}/SGELK";
  for (qw(workingdir numnodes numcpus keep qsubxopts queue)){
    $sge->set($_,$$settings{$_}) if($$settings{$_});
  }
  # Ensure there is no cluster scheduling, if requested
  if(defined($$settings{qsub}) && !$$settings{qsub}){
    $sge->set('scheduler',0);
    $sge->set('qsub',undef);
  }

  # open the log file now that we know the logdir is there
  open($logmsgFh,">$$settings{logdir}/launch_set.log") or die "ERROR: could not open log file $$settings{logdir}/launch_set.log";
  logmsg "======\nOpened logfile at $$settings{logdir}/launch_set.log\n=====";
  logmsg "Called Lyve-SET as follows:\n\n  $invocation\n";

  my $settingsString="";
  $settingsString.=join("\t=\t",$_,$$settings{$_})."\n" for (sort {$a cmp $b} keys(%$settings));
  logmsg "Raw settings are as follows\n$settingsString";

  my $startTimestamp=time();
  logmsg "Lyve-SET started at ".strftime("\%F \%T",localtime());


  #####################################
  # Go through the major steps of SET #
  #####################################
  simulateReads($settings);
  maskReference($ref,$settings);
  indexReference($ref,$settings);
  mapReads($ref,$$settings{readsdir},$$settings{bamdir},$settings);
  variantCalls($ref,$$settings{bamdir},$$settings{vcfdir},$settings);
  annotateVariants($$settings{richseq},$$settings{vcfdir},$settings);
  compareTaxa($ref,$project,$settings); # SNP matrix, alignment, trees
  #####################################
  # Finished major steps              #
  #####################################
  
  # Make an output directory
  # TODO: put 'out' into set_manage.pl
  mkdir "$project/out";
  # Find the output files
  #my $absDir=rel2abs($$settings{msadir}); # save the abs. path
  my $outPrefix="$project/out"; # consistent output naming
  symlink("../msa/out.RAxML_bipartitions","$outPrefix/RAxML.dnd");
  symlink("../msa/out.filteredMatrix.tsv","$outPrefix/snpmatrix.tsv");
  symlink("../msa/out.informative.fasta","$outPrefix/snp.aln.fasta");
  symlink("../msa/out.pairwise.tsv","$outPrefix/pairwise.tsv");
  symlink("../msa/out.pairwiseMatrix.tsv","$outPrefix/pairwise.matrix.tsv");

  my $stopTimestamp=time();
  logmsg "Finished at ".strftime("\%F \%T",localtime());
  my $duration=$stopTimestamp-$startTimestamp;
  my $minutes=int($duration/60);
  my $seconds=$duration % 60;
  logmsg "Duration: $minutes minutes, $seconds seconds";

  return 0;
}

# Pick information about SET from the Makefile where it is stored.
# If overwrite is set to true, then existing settings will be 
# overwritten. Otherwise, settings will be applied if no other
# value was there beforehand.
sub readGlobalSettings{
  my($overwrite,$settings)=@_;
  $settings={} if(!defined($settings));
  # make a copy of settings to avoid ridiculousness
  my %settingsCopy=%$settings;

  my $confPath="$FindBin::RealBin/../config/LyveSET.conf";
  my $cfg=new Config::Simple($confPath) or die "ERROR: could not open $confPath. You might need to run 'make' still.";
  # Find the configuration that was specified
  my $block=$cfg->block("global");
  die "ERROR: configuration [global] was not found in $confPath!" if(!keys(%$block));
  # Put the presets into the current settings
  for(keys(%$block)){
    if($overwrite){
      $settingsCopy{$_}=$$block{$_};
    } else {
      $settingsCopy{$_}=$$block{$_} if(!$settingsCopy{$_});
    }
  }

  return \%settingsCopy;
}

sub richseqToFasta{
  my($infile,$outfile,$settings)=@_;
  logmsg "Converting $infile to $outfile";
  my $seqin=Bio::SeqIO->new(-file=>"$infile");
  my $seqout=Bio::SeqIO->new(-file=>">$outfile");
  while(my $seq=$seqin->next_seq){
    $seq->desc(" ");
    $seqout->write_seq($seq);
  }
  $seqin->close;
  $seqout->close;
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

    # How many simulated reads should there be with a reasonable coverage?
    my $coverage=100;
    my $readLength=250;
    my $approximateGenomeSize = (-s $asm); # estimate based on file size
    if($approximateGenomeSize < $readLength){
      logmsg "WARNING: I am not simulating reads from $asm because it is too small";
      next;
    }
    my $numReads=$approximateGenomeSize * 100 / ($readLength*2);

    # Four successive jobs that depend on each other:
    #   1. Simulate into paired ends
    #   2. Shuffle
    #   3. Gzip
    #   4. Remove simulated unshuffled reads (depends on the shuffle step and not gzip)
    $sge->pleaseExecute("$exec -d 400 -N $numReads -1 $readLength -2 $readLength -e 0.0 -r 0.0 -R 0.0000 -h '$asm' $out1 $out2 && run_assembly_shuffleReads.pl $out1 $out2 > $outFastq && gzip -v '$outFastq' 2>&1 && rm -v $out1 $out2", {jobname=>"simulate$b",numcpus=>1});
  }
  $sge->wrapItUp();
}

sub maskReference{
  my($ref,$settings)=@_;

  my $unmaskedRegions="$$settings{refdir}/unmaskedRegions.bed";
  return $unmaskedRegions if(-e $unmaskedRegions);

  # Make the directory where these masking coordinates go
  my $maskDir="$$settings{refdir}/maskedRegions";
  mkdir($maskDir) if(!-d $maskDir);
  logmsg "Finding regions to mask, and recording them in $maskDir/*";

  # Mark phage locations if they haven't already been marked
  my $phageRegions="$maskDir/phages.bed";
  if($$settings{'mask-phages'} && !-e $phageRegions){
    logmsg "Masking the reference genome for phages";
    $sge->pleaseExecute("set_findPhages.pl --numcpus $$settings{numcpus} $ref > $phageRegions.tmp && mv -v $phageRegions.tmp $phageRegions",{numcpus=>$$settings{numcpus},jobname=>"set_findPhages"});
  }

  $sge->wrapItUp();

  logmsg "DONE with finding regions to mask";
  # END: finding regions to mask
  ##############################

  #########################
  # Find out which regions should be kept in the final pileups
  # and mark it as such in a final file.

  logmsg "Creating a unique list of coordinates to mask";
  # Find out how long each sequence is first.
  my $in=Bio::SeqIO->new(-file=>$ref);
  my %seqLength;
  while(my $seq=$in->next_seq){
    $seqLength{$seq->id}=$seq->length;
  }

  # Find out what the masked regions are in each bed file
  # and find a union of all masked regions, ie, a unique set of points.
  my %maskedRange  = ();
  for my $bed(glob("$maskDir/*.bed")){
    open(BED,$bed) or die "ERROR: could not open $bed for reading: $!";
    while(<BED>){
      chomp;
      my($contig,$start,$stop)=split /\t/;
      $maskedRange{$contig}||=Array::IntSpan->new;
      $maskedRange{$contig}->set_range($start,$stop,1);
    }
    close BED;
  }

  # Write the masked ranges to a file
  my @masked; 
  my $maskedRegions="$$settings{refdir}/maskedRegions.bed";
  open(MASKEDBED,">$maskedRegions") or die "ERROR: could not open $maskedRegions: $!";
  while(my($contig,$RangeObj)=each(%maskedRange)){
    $RangeObj->consolidate(); # merge ranges
    for my $rangeArr($RangeObj->get_range_list()){
      my($min,$max)=@$rangeArr;
      logmsg "Masking ".join("\t",$contig,$min,$max);
      print MASKEDBED join("\t",$contig,$min,$max)."\n";
    }
  }
  close MASKEDBED;

  # Invert the coordinates to make the unmasked coordinates
  my %unmaskedRange=();
  while(my($contig,$RangeObj)=each(%maskedRange)){
    # Make a range spanning the whole contig
    $unmaskedRange{$contig}||=Array::IntSpan->new;
    $unmaskedRange{$contig}->set_range(1,$seqLength{$contig},1);

    # Find where this contig has been masked
    my @maskedRanges=$maskedRange{$contig}->get_range_list();
    my $maskedRanges=$maskedRange{$contig}->get_range_list();
    logmsg "Deleting from $contig: $maskedRanges";
    # Wherever it is masked, remove it from the unmasked objects
    for my $rangeArr(@maskedRanges){
      $unmaskedRange{$contig}->set_range($$rangeArr[0],$$rangeArr[1],undef);
    }
    # Simplify the set of ranges that are unmasked
    $maskedRange{$contig}->consolidate();
  }

  # Write inverted (unmasked) regions to a file
  my $unmaskedRegionsTmp="$unmaskedRegions.tmp";
  open(UNMASKEDBED,">$unmaskedRegionsTmp") or die "ERROR: could not open $unmaskedRegionsTmp: $!";
  while(my($contig,$RangeObj)=each(%unmaskedRange)){
    for my $rangeArr($RangeObj->get_range_list){
      my($min,$max)=@$rangeArr;
      logmsg "Unmasked: ".join("\t",$contig,$min,$max);
      print UNMASKEDBED join("\t",$contig,$min,$max)."\n";
    }
  }
  close UNMASKEDBED;
  system("mv -v $unmaskedRegions.tmp $unmaskedRegions");

  logmsg "Inverted masked regions from $maskedRegions and put them into $unmaskedRegions";
  logmsg "I now know where unmasked regions are, as they are listed in $unmaskedRegions.";

  return $unmaskedRegions;
}

sub indexReference{
  my($ref,$settings)=@_;
  logmsg "Indexing the reference for read mapping";

  # sanity check: see if the reference has reserved characters
  my $in=Bio::SeqIO->new(-file=>$ref);
  while(my $seq=$in->next_seq){
    my $defline=$seq->id." ".$seq->desc;
    die "Offending character found in the defline\n $defline\n $1" if($defline=~/([\:])/);
  }

  logmsg "Indexing with $$settings{mapper}";
  if($$settings{mapper} eq 'smalt'){
    return $ref if(-e "$ref.sma" && -e "$ref.smi");
    system("smalt index -k 13 -s 2 $ref $ref 2>&1");
    die if $?;
  } elsif($$settings{mapper} eq 'snap'){
    return $ref if(-d "$ref.snap" && -e "$ref.snap/GenomeIndex");
    # -bSpace to make sure the deflines match up with what samtools expects
    # -exact -large to take the time to make a fast but somewhat smaller db
    # Must specify the number of CPUs or else it greedily takes them all instead of what the user specifies
    system("snap index $ref $ref.snap -exact -large -bSpace -t$$settings{numcpus}");
    die if $?;
  } elsif($$settings{mapper} eq 'stampy'){
    return $ref if(-e "$ref.stidx" && -e "$ref.sthash");
    system("stampy.py --noparseNCBI -G $ref $ref && stampy.py --noparseNCBI -g $ref -H $ref");
    die if $?;
  }
  elsif($$settings{mapper} eq 'bowtie2'){
    return $ref if(-e "$ref.{1,2,3,4}.bt2" && -e "$ref.rev.{1,2}.bt2");
    system("bowtie2-build -q -f $ref $ref");
    die if $?;
  }
  elsif($$settings{mapper} eq 'bwa'){
    return $ref if(-e "$ref.{amb,ann,bwt,pac,sa}");
    system("bwa index $ref");
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
  # TODO use LyveSET.pm fastq extensions for this step, to streamline
  my @file=(glob("$readsdir/*.fastq"),glob("$readsdir/*.fastq.gz"),glob("$readsdir/*.fq"),glob("$readsdir/*.fq.gz"));

  if(!@file){
    logmsg "ERROR: no files were found in $readsdir/. However, in case you are continuing Lyve-SET after all the bam files have already been created and you just happened to delete all the input files, I will skip read mapping and will continue";
    return 0;
  }

  downsampleReads($ref,\@file,$settings) if($$settings{downsample});
  cleanReads(\@file,$settings) if($$settings{read_cleaner});

  my $snapxopts="";
  my $smaltxopts="";
  my $stampyxopts="";
  my $bowtie2xopts="";
  my $bwaxopts="";
  if($$settings{singleend}){
    $smaltxopts.="--pairedend 1";
    $snapxopts.= "--pairedend 1";
    $stampyxopts.= "--pairedend 1";
    $bwaxopts.= "--pairedend 1";
  }

  my (@bam,@job);
  for my $fastq(@file){
    my $b=fileparse $fastq;
    my $bamPrefix="$bamdir/$b-".basename($ref,@fastaExt);
    push(@bam,"$bamPrefix.sorted.bam");

    if(-e "$bamPrefix.sorted.bam"){
      logmsg "Found $bamPrefix.sorted.bam. Skipping.";
      next;
    }else{
      logmsg "Mapping to create $bamPrefix.sorted.bam";
    }

    if($$settings{mapper} eq 'smalt'){
      $sge->pleaseExecute("$scriptsdir/launch_smalt.pl $smaltxopts -ref $ref -f $fastq -b $bamPrefix.sorted.bam -tempdir $tmpdir --numcpus $$settings{numcpus} ",{jobname=>"smalt$b"});
    } elsif($$settings{mapper} eq 'snap'){
      $sge->pleaseExecute("$scriptsdir/launch_snap.pl $snapxopts -ref $ref -f $fastq -b $bamPrefix.sorted.bam -tempdir $tmpdir --numcpus $$settings{numcpus} ",{jobname=>"snap$b"});
    } elsif($$settings{mapper} eq 'stampy'){
      # Make a tmpdir for stampy since each invocation needs its own space
      my $stampydir=tempdir("$tmpdir/stampyXXXXXX",CLEANUP=>1);
      $sge->pleaseExecute("$scriptsdir/launch_stampy.sh $stampyxopts -r $ref -f $fastq -b $bamPrefix.sorted.bam -t $stampydir --numcpus $$settings{numcpus} ",{jobname=>"stampy$b"});
    } 
    elsif($$settings{mapper} eq 'bowtie2'){
      my $bowtie2dir=tempdir("$tmpdir/bowtie2XXXXXX",CLEANUP=>1);
      $sge->pleaseExecute("$scriptsdir/launch_bowtie2.sh $bowtie2xopts -r $ref -f $fastq -b $bamPrefix.sorted.bam -t $bowtie2dir --numcpus $$settings{numcpus} ",{jobname=>"bowtie2$b"});
    }
    elsif($$settings{mapper} eq 'bwa'){
      my $bwadir=tempdir("$tmpdir/bwaXXXXXX",CLEANUP=>1);
      $sge->pleaseExecute("$scriptsdir/launch_bwa.sh $bwaxopts -r $ref -f $fastq -b $bamPrefix.sorted.bam -t $bwadir --numcpus $$settings{numcpus} ",{jobname=>"bwa$b"});
    }
    else {
      die "ERROR: I do not understand the mapper $$settings{mapper}";
    }
  }
  logmsg "All mapping jobs have been submitted. Waiting on them to finish.";
  $sge->wrapItUp();

  # Look for cliffs in coverage levels
  if($$settings{'mask-cliffs'}){
    # Make the directory where these masking coordinates go
    my $maskDir="$$settings{refdir}/maskedRegions";
    mkdir($maskDir) if(!-d $maskDir);
    logmsg "Finding regions to mask, and recording them in $maskDir/*";

    for my $bam(@bam){
      my $b=basename($bam,@bamExt);
      my $bed="$bam.cliffs.bed";
      next if(-e $bed);
      $sge->pleaseExecute("$scriptsdir/set_findCliffs.pl --numcpus $$settings{numcpus} $bam > $bed",{jobname=>"maskCliff$b",numcpus=>$$settings{numcpus}});
    }
    $sge->wrapItUp();
  }

  return 1;
}

sub downsampleReads{
  my($ref,$reads,$settings)=@_;
  my $coverage=50;
  my $reducedBases=(-s $ref) * $coverage;
  logmsg "Downsampling the reads to ${coverage}x.";
  for my $file(@$reads){
    # downsample into tmpdir
    # move the original reads to "$file.orig"
    my $b=basename($file,qw(.cleaned.fastq.gz .downsampled.fastq.gz),@fastqExt);
    my $d=dirname($file);
    my $backupDir="$d/notDownsampled";
    mkdir($backupDir);
    my $backup="$backupDir/$b.fastq.gz";
    next if(-e $backup);
    logmsg "Did not find $backup. Downsampling...";
    # 1. Downsample to file.downsampled
    # 2. Mv un-downsampled file to file.orig
    # 3. Mv file.downsampled to file
    $sge->pleaseExecute("run_assembly_removeDuplicateReads.pl -size $reducedBases $file | gzip -c > $file.downsampled && mv -v $file $backup && mv -v $file.downsampled $file",{numcpus=>1,jobname=>"downsample$b"});
  }
  $sge->wrapItUp();
}

sub cleanReads{
  my($reads,$settings)=@_;
  for my $file(@$reads){
    # Clean into tmpdir
    # move the original reads to "$file.orig"
    my $b=basename($file,qw(.cleaned.fastq.gz .downsampled.fastq.gz), @fastqExt);
    my $d=dirname($file);
    my $backupDir="$d/notCleaned";
    mkdir($backupDir);
    my $backup="$backupDir/$b.fastq.gz";
    my $tmp="$$settings{tmpdir}/$b.fastq.gz";
    next if(-e $backup);
    # 1. Clean to tmp/cleaned
    # 2. Mv uncleaned to backup/file.fastq.gz
    # 3. mv tmp/cleaned to ./file.fastq.gz
    if($$settings{read_cleaner} eq "cgp"){
      logmsg "Did not find $backup. Cleaning with CGP...";
      $sge->pleaseExecute("run_assembly_trimClean.pl -i $file -o $tmp --nosingletons --numcpus $$settings{numcpus} $$settings{read_cleaner_cgpxopts} 2>&1 && mv -v $file $backup && mv -v $tmp $file",{numcpus=>$$settings{numcpus},jobname=>"trimClean$b"});
    } elsif($$settings{read_cleaner} eq "bayeshammer"){
      logmsg "Did not find $backup. Cleaning with BayesHammer...";
      my $tmpdir=tempdir("$$settings{tmpdir}/bayeshammerXXXXXX",CLEANUP=>1);
      $sge->pleaseExecute("set_bayesHammer.pl $file --numcpus $$settings{numcpus} --tempdir $tmpdir | gzip -c > $tmp && mv -v $file $backup && mv -v $tmp $file",{numcpus=>$$settings{numcpus},jobname=>"bayeshammer$b"});
    } else {
      logmsg "ERROR: I do not understand the read cleaner $$settings{read_cleaner}";
      die;
    }
  }
  $sge->wrapItUp();
}

sub variantCalls{
  my($ref,$bamdir,$vcfdir,$settings)=@_;
  logmsg "Calling variants with $$settings{snpcaller}";
  my @bam=glob("$bamdir/*.sorted.bam");

  if(!@bam){
    logmsg "ERROR: no bam files were found in $bamdir/. However, in case you are continuing Lyve-SET after all the vcf files have already been created and you just happened to delete all the bam files, I will continue";
    return 0;
  }

  # see if bgzip and tabix exist
  my $bgzip=`which bgzip 2>/dev/null`; chomp($bgzip);
  my $tabix=`which tabix 2>/dev/null`; chomp($tabix);

  # If --fast is given, choose a genome. Find possible 
  # variant sites in that genome and use those regions 
  # for the next SNP-calling instances.
  my $regionsFile="";
  if($$settings{'sample-sites'}){
    logmsg "--sample-sites was specified: finding initial variant sites to make the downstream SNP calling faster. First, need to find most distant genome using kmers and jaccard distance";

    # Find the most different set of reads.
    # TODO: just find the max distance from the ref assembly.
    for(my $i=0;$i<@bam;$i++){
      for(my $j=$i+1;$j<@bam;$j++){
        my $tmpName="$$settings{tmpdir}/".join("_","genomeDist",$i,$j,".tmp");
        $sge->pleaseExecute("$scriptsdir/genomeDist.pl -n 1 $bam[$i] $bam[$j] > $tmpName",{numcpus=>1,jobname=>"genomeDist"});
      }
    }
    $sge->wrapItUp();

    # Figure out the most distant using sums of points.
    my %distance;
    open(SORTALLGENOMEDIST,"sort -k3,3n -k1,1 -k2,2 $$settings{tmpdir}/genomeDist*.tmp | ") or die "ERROR: could not open  $$settings{tmpdir}/genomeDist*.tmp: $!";
    while(<SORTALLGENOMEDIST>){
      chomp;
      my($bam1,$bam2,$d)=split /\t/;
      $distance{$bam1}+=$d;
      $distance{$bam2}+=$d;
    }
    close SORTALLGENOMEDIST;

    # The most distant genome is the one with the lowest jaccard distance overall.
    my $initBam=$bam[0];
    my $minScore=$distance{$initBam};
    for my $b(@bam){
      if($distance{$b} > $minScore){
        $minScore=$distance{$b};
        $initBam=$b;
      }
    }
    logmsg "The most distantly related by jaccard distance is $initBam.";
    system("rm -vf $$settings{tmpdir}/genomeDist*.tmp");

    # This complicated command does multithreaded samtools mpileup
    # using xargs and one cpu per contig.
    # https://www.biostars.org/p/48781/#48815
    $regionsFile="$$settings{tmpdir}/initialSnpSites.bed";
    if(-e $regionsFile && !-s $regionsFile){
      unlink $_ for(glob("$regionsFile*"));
    }
    if(!-e $regionsFile){
      system("rm -fv $regionsFile.*.tmp");
      my $command="samtools view -H $initBam | grep \"\@SQ\" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P $$settings{numcpus} sh -c \"samtools mpileup -suf $ref -r '{}' $initBam | bcftools call -c | grep -v '^#' | cut -f 1,2 > $regionsFile.'{}'.tmp \" ";
      $sge->pleaseExecute($command,{numcpus=>$$settings{numcpus},jobname=>"snpPositionDiscovery"});
      $sge->wrapItUp();

      system("sort -k1,1 -k2,2n $regionsFile.*.tmp | uniq | grep . > $regionsFile && rm -v $regionsFile.*.tmp");
      die if $?;
    }
  }
  # # In at least one rare event, sampling sites introduced a zero-byte file
  # # which made Lyve-SET crash.  Avoid, avoid, avoid!
  # if(-e $regionsFile && !-s $regionsFile){
  #   unlink $_ for(glob("$regionsFile*"));
  # }

  # TODO make the variant callers output to bgzip and indexed files
  for my $bam(@bam){
    my $b=fileparse($bam,".sorted.bam");
    if(-e "$vcfdir/$b.vcf.gz.tbi"){
      logmsg "Found $vcfdir/$b.vcf.gz. Skipping";
      next;
    }
    logmsg "Calling SNPs into $vcfdir/$b.vcf";
    my $jobname=""; # the job that the filter/sort scripts are waiting for
    if($$settings{snpcaller} eq 'varscan'){
      $jobname="varscan$b";
      my $varscanxopts="";
      $varscanxopts.="--region $regionsFile " if($regionsFile);
      $varscanxopts.="--exclude $bam.cliffs.bed " if(-s "$bam.cliffs.bed");
      my $varscanCommand="$scriptsdir/launch_varscan.pl $bam --numcpus $$settings{numcpus} --tempdir $$settings{tmpdir} --reference $ref --altfreq $$settings{min_alt_frac} --coverage $$settings{min_coverage} $varscanxopts > $vcfdir/$b.vcf";
      logmsg $varscanCommand;
      $sge->pleaseExecute($varscanCommand,{numcpus=>$$settings{numcpus},jobname=>$jobname,qsubxopts=>""});
      # sort VCF
      $sge->pleaseExecute("mv $vcfdir/$b.vcf $vcfdir/$b.vcf.tmp && 
          vcf-sort < $vcfdir/$b.vcf.tmp > $vcfdir/$b.vcf && 
          rm -v $vcfdir/$b.vcf.tmp",
          {jobname=>"sort$b",qsubxopts=>"-hold_jid $jobname",numcpus=>1}
      );
      $jobname="sort$b"; # the job that bgzip waits on to finish

    }elsif($$settings{snpcaller} eq 'vcftools'){
      $jobname="vcftools$b";
      my $vcftoolsxopts="";
      $vcftoolsxopts.="--region $$settings{refdir}/unmaskedRegions.bed" if(-s "$$settings{refdir}/unmaskedRegions.bed");
      $vcftoolsxopts.="--exclude $bam.cliffs.bed " if(-s "$bam.cliffs.bed");
      my $vcftoolsCommand="$scriptsdir/launch_vcftools.sh -m $bam --numcpus $$settings{numcpus} --tempdir $$settings{tmpdir} --reference $ref --altfreq $$settings{min_alt_frac} --coverage $$settings{min_coverage} $vcftoolsxopts";
      logmsg $vcftoolsCommand;
      $sge->pleaseExecute($vcftoolsCommand,{numcpus=>$$settings{numcpus},jobname=>$jobname,qsubxopts=>""});
      $sge->pleaseExecute("mv -fv $$settings{tmpdir}/$b.vcf.gz $vcfdir/$b.vcf.gz &&
        mv -fv $$settings{tmpdir}/$b.vcf.gz.tbi $vcfdir/$b.vcf.gz.tbi",
        {jobname=>"moveVCF$b",qsubxopts=>"-hold_jid $jobname",numcpus=>1}
        );
      $jobname="moveVCF$b"; # the job that waits on to finish
    }

    else {
      die "ERROR: I do not understand snpcaller $$settings{snpcaller}";
    }

    # bgzip and tabix indexing
    indexAndCompressVcf("$vcfdir/$b.vcf",$jobname,$settings);

  } # END bam while loop

  logmsg "All variant-calling jobs have been submitted. Waiting on them to finish";
  $sge->wrapItUp();

  return 1;
}


sub annotateVariants{
  my($richseq,$vcfdir,$settings)=@_;
  return if(!$richseq);

  for my $vcf(glob("$vcfdir/*.vcf.gz")){
    my $b=basename($vcf,@vcfExt);
    $sge->pleaseExecute("launch_snpeff.pl --genbank $richseq --outdir $vcfdir.annotated --numcpus 1 $vcf",{numcpus=>1,jobname=>"snpeff$b"});
  }
  # No need to wrapItUp() since nothing depends on these annotations
}


# returns the last SGE job information in a hash ref
sub indexAndCompressVcf{
  my($vcf,$holdjid,$settings)=@_;
  my $j={};
  eval{
    $j=$sge->pleaseExecute("
      vcf-sort < '$vcf' > '$vcf.sorted.tmp' && mv '$vcf.sorted.tmp'  '$vcf' && \
      bgzip -f '$vcf' && tabix '$vcf.gz'
    ",{qsubxopts=>"-hold_jid $holdjid",jobname=>"sortAndCompress",numcpus=>1});
  };
  if($@){
    logmsg "Warning: Could not compress and index $vcf: $@\n  See if bgzip and tabix are both installed correctly.";
    die;
  }
  return $j;
}

# Make a matrix, trees, etc
sub compareTaxa{
  my($ref,$project,$settings)=@_;
  my $pooled;
  if($$settings{matrix}){
    $pooled=variantsToMatrix($ref,$$settings{bamdir},$$settings{vcfdir},$$settings{msadir},$settings);
  } else {
    logmsg "The matrix was not requested; wrapping up";
    return 1;
  }

  if($$settings{msa} || $$settings{trees}){
    logmsg "Launching set_processPooledVcf.pl";
    my $command="set_processPooledVcf.pl $pooled --allowedFlanking $$settings{allowedFlanking} --prefix $$settings{msadir}/out --numcpus $$settings{numcpus} --exclude $$settings{refdir}/maskedRegions.bed 2>&1 | tee --append $$settings{logdir}/launch_set.log";
    logmsg "Processing the pooled VCF\n  $command";
    $sge->pleaseExecute($command,{numcpus=>$$settings{numcpus},jobname=>"set_processPooledVcf.pl"});
    $sge->wrapItUp();
  } else {
    logmsg "The phylogeny was not requested; wrapping up";
  }

  # Diagnose the project
  # Make sure this gets printed to stdout
  logmsg "Running set_diagnose.pl";
  #my $diagnosis="$$settings{logdir}/diagnosis.txt";
  $sge->pleaseExecute("set_diagnose.pl -p $project 2>&1 | tee --append $$settings{logdir}/launch_set.log",{jobname=>"set_diagnose",numcpus=>1});
  $sge->wrapItUp();

  return 1;
}


sub variantsToMatrix{
  my ($ref,$bamdir,$vcfdir,$msadir,$settings)=@_;
  logmsg "Creating a core hqSNP matrix";
  my $logdir=$$settings{logdir};

  my $matrix="$msadir/out.matrix.tsv";
  my $unFixedPooled="$msadir/out.unfixed.pooled.vcf.gz";
  my $unFixedSnpPooled="$msadir/out.unfixed.pooled.snps.vcf.gz";
  my $pooled="$msadir/out.pooled.vcf.gz";
  my $snpPooled="$msadir/out.pooled.snps.vcf.gz";
  if(-e $pooled){
    logmsg "Found $pooled -- already present. Not re-converting.";
    return $pooled;
  }

  # input files
  my @vcf=glob("$vcfdir/*.vcf.gz");
  my @bam=glob("$bamdir/*.sorted.bam");
  my @infile=(@bam,@vcf);
  my $infile="'".join("' '",@infile)."'";
  my $inVcf="'".join("' '",@vcf)."'";

  # mergevcf deletes its tmpdir and so we need to make a new one
  my $tmpdir=tempdir("$$settings{tmpdir}/mergevcfXXXXXX",CLEANUP=>0);

  my $mergexopts="";
  $mergexopts.="-r $ref.regions.txt" if(-e "$ref.regions.txt" && -s "$ref.regions.txt" > 0);
  my $mergeCommand="mergeVcf.sh $mergexopts -s -n $$settings{numcpus} -t $tmpdir -o $unFixedPooled $inVcf >&2";
  logmsg $mergeCommand;
  $sge->pleaseExecute($mergeCommand,{jobname=>"poolVcfs",numcpus=>$$settings{numcpus}});
  $sge->wrapItUp();

  # Spruce up the VCF so that it conforms to Lyve-SET thresholds.
  logmsg "Reevaluating the pooled VCF with set_fixVcf.pl";
  $sge->pleaseExecute("set_fixVcf.pl --fail-samples --pass-until-fail --DP4 2 --min_coverage $$settings{min_coverage} --min_alt_frac $$settings{min_alt_frac} $unFixedPooled > $$settings{tmpdir}/out.pooled.vcf && bgzip -c $$settings{tmpdir}/out.pooled.vcf > $pooled && tabix -f $pooled",{jobname=>"fixPooled",numcpus=>1});
  $sge->pleaseExecute("set_fixVcf.pl --fail-samples --pass-until-fail --DP4 2 --min_coverage $$settings{min_coverage} --min_alt_frac $$settings{min_alt_frac} $unFixedSnpPooled> $$settings{tmpdir}/out.snps.vcf && bgzip -c $$settings{tmpdir}/out.snps.vcf > $snpPooled && tabix -f $snpPooled",{jobname=>"fixSnpsPooled",numcpus=>1});
  $sge->wrapItUp();

  logmsg "Done reevaluating.";
  system("rm -f $unFixedPooled $unFixedSnpPooled");

  return $pooled;
}

sub usage{
  my($settings)=@_;

  ## Format a few variables correctly
  # simplified pathnames for some options
  my @dir=qw(asmdir msadir readsdir bamdir vcfdir tmpdir logdir);
  $$settings{$_}=abs2rel($_).'/' for(@dir);
  # right padding for some options
  $$settings{$_}=reverse(sprintf("%15s","".reverse($$settings{$_}))) for(qw(mapper snpcaller allowedFlanking),@dir);
  local $0=fileparse $0;

  # The help menu
  my $help="Launches the Lyve-SET pipeline
    Please visit http://github.com/lskatz/Lyve-SET for more information.

    Usage: $0 [project] [-ref reference.fasta|reference.gbk]

    If project is not given, then it is assumed to be the current working directory.
    If reference is not given, then it is assumed to be proj/reference/reference.fasta
    -ref      proj/reference/reference.fasta   The reference genome assembly. If it is 
                                               a genbank or embl file, then it will be 
                                               converted to reference.gbk.fasta and will
                                               be used for SNP annotation. If a fasta
                                               is given, then no SNP annotation will
                                               happen. Using a gbk or embl file is currently
                                               experimental.

    SNP FILTER OPTIONS
    --allowedFlanking  $$settings{allowedFlanking}allowed flanking distance in bp. Nucleotides this close together cannot be considered as high-quality.
    --min_alt_frac     $$settings{min_alt_frac}        The percent consensus that needs to be reached before a SNP is called. Otherwise, 'N'
    --min_coverage     $$settings{min_coverage}          Minimum coverage needed before a SNP is called. Otherwise, 'N'
    ";
    return "$help\n  --help To view more help\n\n" if(!$$settings{help});

    $help.="
    Where parameters with a / are directories
    LOCATIONS OF FILE DIRECTORIES
    -reads      $$settings{readsdir}    where fastq and fastq.gz files are located
    --bamdir    $$settings{bamdir}    where to put bams
    --vcfdir    $$settings{vcfdir}    where to put vcfs
    --tmpdir    $$settings{tmpdir}    tmp/ Where to put temporary files
    --msadir    $$settings{msadir}    multiple sequence alignment and tree files (final output)
    --logdir    $$settings{logdir}    Where to put log files. Qsub commands are also stored here.
    --asmdir    $$settings{asmdir}    directory of assemblies. Copy or symlink the reference genome assembly 
                                   to use it if it is not already in the raw reads directory

    PERFORM CERTAIN STEPS
    --mask-phages                  Search for and mask phages in the reference genome
    --mask-cliffs                  Search for and mask 'Cliffs' in pileups

    SKIP CERTAIN STEPS
    --nomatrix                     Do not create an hqSNP matrix
    --nomsa                        Do not make a multiple sequence alignment
    --notrees                      Do not make phylogenies
    --singleend                    Treat everything like single-end. Useful for 
                                   when you think there is a single-end/paired-end bias.
    OTHER SHORTCUTS
    --fast                         Shorthand for --downsample --mapper snap --nomask-phages 
                                                 --nomask-cliffs --sample-sites
    --presets \"\"                   See presets.conf for more information
    --downsample                   Downsample all reads to 50x. Approximated according 
                                   to the ref genome assembly
    --sample-sites                 Randomly choose a genome and find SNPs in a quick 
                                   and dirty way. Then on the SNP-calling stage, 
                                   only interrogate those sites for SNPs for each 
                                   genome (including the randomly-sampled genome).
    MODULES
    --read_cleaner none            Which read cleaner? Choices: none, CGP, BayesHammer
    --mapper       $$settings{mapper} Which mapper? Choices: smalt, snap, stampy, bowtie2, bwa
    --snpcaller    $$settings{snpcaller} Which SNP caller? Choices: varscan, vcftools

    SCHEDULER AND MULTITHREADING OPTIONS
    --queue        $$settings{queue}           default queue to use
    --numnodes     $$settings{numnodes}              maximum number of computer nodes on the scheduler
    --numcpus      $$settings{numcpus}               number of cpus per node
    --qsubxopts    '-N lyve-set'   Extra options to pass to qsub. This is not 
                                   sanitized; internal options might overwrite yours.
    --noqsub                       Do not use the scheduler, even if it exists
  ";
  return $help;
}
