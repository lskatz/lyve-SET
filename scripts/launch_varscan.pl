#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Creates a vcf using varscan

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename fileparse/;
use File::Temp qw/tempdir/;
use Bio::Perl;
use Bio::FeatureIO;

use threads;
use threads::shared;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib/vcftools_0.1.12b/perl";
use lib "$FindBin::RealBin/../lib";
use Vcf;
use LyveSET qw/@fastaExt @fastqExt @bamExt logmsg/;

my $samplename="Sample1"; # to be changed later
my $readBamStick :shared; # only one bam can be read at a time
$0=fileparse $0;
exit(main());
 

sub main{
  my $settings={};
  GetOptions($settings,qw(help reference=s tempdir=s altFreq=s coverage=i region=s exclude=s numcpus=i)) or die $!;

  my $bam=$ARGV[0] or die "ERROR: need bam\n".usage();
  die "ERROR: only one bam allowed right now" if (@ARGV > 1);
  my $reference=$$settings{reference} or die "ERROR: need --reference\n".usage();
  my $refname=basename($reference,@fastaExt);

  $$settings{tempdir} ||=tempdir("set_varscan.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{coverage}||=10;
  $$settings{altFreq} ||=0.75;
  $$settings{region}  ||="";
  $$settings{exclude} ||="";
  $$settings{numcpus} ||=1;

  ## Parameter error checking
  # Not sure if I should include this error check or not
  # die "ERROR: Cannot use --region and --exclude together\n".usage() if($$settings{exclude} && $$settings{region});
  die "ERROR: altFreq must be between 0 and 1, inclusive" if($$settings{altFreq} < 0 || $$settings{altFreq} > 1);
  die "ERROR: coverage must be >=0" if($$settings{coverage} < 0);

  # Check to see if varscan is installed correctly
  `varscan.sh >/dev/null`;
  die "ERROR: varscan.sh gave an error" if $?;

  # Transform the sample name, e.g.,
  # NC001416.fasta.wgsim.fastq.gz-reference.sorted.bam => NC001416.fasta.wgsim
  # sample1.fastq.gz-reference.sorted.bam => sample1
  $samplename=basename($bam,@bamExt);                    # remove bam extension
  $samplename=~s/\-$refname$//;                          # remove reference genome name
  $samplename=basename($samplename,@fastqExt,@fastaExt); # remove fast[aq] extension
  $samplename=basename($samplename,@fastaExt);           # remove fasta extension

  my $vcf=varscan2($bam,$reference, $settings);

  # Print the vcf to stdout
  open(VCF,$vcf) or die "ERROR: could not read $vcf: $!";
  print while(<VCF>);
  close VCF;

  # remove temporary files
  #unlink($_) for($pileup,$vcf);

  return 0;
}

sub varscan2{
  my($bam,$reference,$settings)=@_;
  # I like to use the shorthand $b for basename, but it is a reserved
  # variable name in sort{} and so I need to use 'local' instead of 'my'
  local $b=fileparse($bam);
  my $pileup="$$settings{tempdir}/$b.mpileup";
  return $pileup if(-e $pileup && -s $pileup > 0);

  # Make sure this whole operation is self-contained
  # and that the temporary dir exists for this sample.
  my $tmpdir="$$settings{tempdir}/$samplename";
  system("rm -rf $tmpdir");
  mkdir($tmpdir);
  # Make sure that no other pileup file gets in the way.
  system("rm -vf $$settings{tempdir}/$b*.mpileup $$settings{tempdir}/$b*.vcf.gz* >&2");

  # Kick off multithreading
  my $regionQueue=Thread::Queue->new();
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&varscanWorker,$bam,$reference,$regionQueue,$settings);
  }

  # Figure out any regions
  my @regions=`makeRegions.pl --numcpus $$settings{numcpus} --numchunks $$settings{numcpus} $bam`;
  chomp(@regions);
  # Add the regions and some thread-terminators to the queue
  $regionQueue->enqueue(@regions);
  $regionQueue->enqueue(undef) for(@thr);
  # When each thread encounters an undef, it will terminate.
  # When the threads are terminated and joined, we
  # can continue past this line.
  for(@thr){
    $_->join;
    die "ERROR with thread ".$_->tid if($_->error());
  }

  # Sort all subVCFs by their starting position since each one is 
  # already internally sorted.
  my @subVcf=glob("$tmpdir/*.vcf.gz");
  my $subVcf=join(" ",@subVcf);

  # Join all pieces of the sample's vcf. All pieces
  # have been sorted at this point.
  logmsg "Concatenating all the VCFs";
  system("bcftools concat --allow-overlaps --remove-duplicates $subVcf > $$settings{tempdir}/$b.merged.vcf");
  die if $?;

  # Make sure everything is cleaned up whenever the script ends
  system("rm -rf $tmpdir");

  return "$$settings{tempdir}/$b.merged.vcf";
}
  

sub varscanWorker{
  my($bam,$reference,$regionQueue,$settings)=@_;

  my $TID=threads->tid;

  # Figure out mpileup options
  my $mpileupxopts="";
  $mpileupxopts.="-Q 0 -B "; # assume that all reads have been properly filtered at this point and that mappings are good
  $mpileupxopts.="--positions $$settings{region} " if($$settings{region});

  my $tmpdir=tempdir("$$settings{tempdir}/$samplename/varscanWorker.XXXXXX");

  # Make a sample list for varscan so that it doesn't say "Sample1"
  my $sampleList="$tmpdir/samples.txt";
  open(SAMPLELIST,">",$sampleList) or die "ERROR: could not write to $sampleList: $!";
  print SAMPLELIST $samplename;
  close SAMPLELIST;

  my $vcfCounter=0;
  while(defined(my $region=$regionQueue->dequeue)){
    my $pileup="$tmpdir/TID$TID.$vcfCounter.pileup";
    my $vcf="$tmpdir/TID$TID.$vcfCounter.vcf";
    # Avoid disk I/O problems by sleeping
    # a bit of time, depending on the number of threads
    #sleep ($TID % $$settings{numcpus});

    # Avoid disk I/O problems by doing only one mpileup at a time
    {
      lock $readBamStick;

      # Run mpileup on this one region
      my $mpileupCommand="samtools mpileup -f $reference $mpileupxopts --region '$region' $bam > $pileup";
      system($mpileupCommand);
      die "ERROR running mpileup on region $region:\n  $mpileupCommand" if $?;
    }

    # Estimate the size of the pileup. If it's zero, then
    # there will be no SNPs in it and should be skipped.
    if( -s $pileup){
      logmsg "There is a pileup in region $region and so I will run varscan";
    } else {
      logmsg "No pileup was formed in region $region and so there will be no SNPs. Skipping.";
      next;
    }

    # Finally run varscan on the pileup. Bgzip and index the
    # output file. Later it will be combined with bcftools concat.
    my $varscanCommand="varscan.sh mpileup2cns $pileup --min-coverage $$settings{coverage} --min-var-freq $$settings{altFreq} --output-vcf 1 --min-avg-qual 0 --vcf-sample-list $sampleList > $vcf.tmp.vcf";
    system($varscanCommand);
    die "ERROR with varscan!\n  $varscanCommand" if $?;

    # Fix the Vcf
    my $fixedVcf=backfillVcfValues("$vcf.tmp.vcf",$settings);

    my $newVcfLoc="$$settings{tempdir}/$samplename/".basename($vcf);
    system("mv -v $fixedVcf $newVcfLoc >&2"); die if $?;
    system("bgzip $newVcfLoc && tabix $newVcfLoc.gz");
    die if $?;
    $vcfCounter++;
  }

  # Some cleanup
  system("rm -rf $tmpdir");
}

sub mpileup{
  my($bam,$reference,$settings)=@_;
  my $b=fileparse($bam);
  my $pileup="$$settings{tempdir}/$b.mpileup";
  return $pileup if(-e $pileup && -s $pileup > 0);
  logmsg "Creating a pileup $pileup";

  # Figure out any defined seqnames
  my @seqname;
  open(SAMTOOLS,"samtools view -H '$bam' | ") or die "ERROR: could not open $bam with samtools: $!";
  while(<SAMTOOLS>){
    next if(!/^\@SQ/);
    chomp;
    my @F=split /\t/;
    for(@F){
      my($key,$value)=split /:/;
      push(@seqname,$value) if($key eq 'SN');
    }
  }
  close SAMTOOLS;
  my $regions=join("\n",@seqname)."\n";

  # Figure out mpileup options
  my $xopts="";
  $xopts.="-Q 0 -B "; # assume that all reads have been properly filtered at this point and that mappings are good
  $xopts.="--positions $$settings{region} " if($$settings{region});

  # Multithread a pileup so that each region gets piped into a file.
  # Then, these individual files can be combined at a later step.
  system("rm -fv $$settings{tempdir}/$b.*.mpileup >&2"); # get rid of any files that might be in the way
  my $command="echo \"$regions\" | xargs -P $$settings{numcpus} -n 1 -I $0 sh -c 'echo \"MPileup on $0\" >&2; samtools mpileup -f $reference $xopts --region \"$0\" $bam > $$settings{tempdir}/$b.\$\$.mpileup' ";
  logmsg "Running mpileup:\n  $command";
  system($command);
  die "ERROR with xargs and samtools mpileup" if $?;

  # Concatenate the output into a single file.
  # Individual regions are already sorted, thanks to the way mpileup works.
  logmsg "Sorting mpileup results into a combined file";
  system("cat $$settings{tempdir}/$b.*.mpileup > $pileup && rm -f $$settings{tempdir}/$b.*.mpileup");
  die "ERROR with sorting mpileup results and then deleting the intermediate files" if $?;

  return $pileup;
}

sub varscan{
  my($pileup,$settings)=@_;
  die "ERROR: the pileup is a zero-byte file\n  $pileup" if(-s $pileup < 1);
  my $vcf="$$settings{tempdir}/".fileparse($pileup).".tmp.vcf.gz";
  system("varscan.sh mpileup2cns $pileup --min-coverage $$settings{coverage} --min-var-freq $$settings{altFreq} --output-vcf 1 --min-avg-qual 0 |\
    perl -lane 's/Sample1/\Q$samplename\E/; print;' |\
    bgzip -c > $vcf
  ");
  die "ERROR: problem with either varscan.sh, bgzip, or perl one-liner for regex" if $?;

  return $vcf;
}

sub backfillVcfValues{
  my($vcfFile,$settings)=@_;

  my $excludeSites=readBed($$settings{exclude},$settings) if($$settings{exclude});

  my $newVcf="$vcfFile.backfilled.vcf";

  logmsg "Backfilling values in $vcfFile and printing to $newVcf";
  open(my $vcfOut, ">", $newVcf) or die "ERROR: could not write to $newVcf";
  my $vcf=Vcf->new(file=>$vcfFile);

  # Reasons a SNP doesn't pass
  my $fail_lowCoverage="DP$$settings{coverage}";
  my $fail_lowRefFreq ="RF$$settings{altFreq}";
  my $fail_lowAltFreq ="AF$$settings{altFreq}";
  my $fail_indel      ="isIndel";
  my $fail_masked     ="masked";

  # Add new headers
  $vcf->add_header_line({key=>'reference',value=>$$settings{reference}});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowCoverage, Description=>"Depth is less than $$settings{coverage}, the user-set coverage threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowRefFreq, Description=>"Reference variant consensus is less than $$settings{altFreq}, the user-set threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowAltFreq, Description=>"Allele variant consensus is less than $$settings{altFreq}, the user-set threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_indel, Description=>"Indels are not used for analysis in Lyve-SET"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_masked, Description=>"This site was masked using a bed file or other means"});

  # Genotype fields
  $vcf->add_header_line({key=>'FORMAT', ID=>'AF', Number=>'A', Type=>'Float', Description=>"Allele Frequency"});
  $vcf->add_header_line({key=>'FORMAT', ID=>'RF', Number=>'A', Type=>'Float', Description=>"Reference Frequency"});
  $vcf->add_header_line({key=>'FORMAT', ID=>'ADP', Number=>'A', Type=>'Integer', Description=>"Depth of bases with Phred score >= 15"});
  $vcf->add_header_line({key=>'FORMAT', ID=>'FT', Number=>1, Type=>'String', Description=>"Genotype filters using the same codes as the FILTER data element"});

  # Remove unwanted headers
  # ADP was removed as an info field but added as a format field
  $vcf->remove_header_line(key=>'INFO',ID=>'ADP');

  # Done with headers; parse them
  $vcf->parse_header();

  # assume there is only one sample name
  my $samplename=($vcf->get_samples)[0];

  # start printing
  print $vcfOut $vcf->format_header();
  while(my $x=$vcf->next_data_hash()){
    my $posId=$$x{CHROM}.':'.$$x{POS};
    my $pos=$$x{POS};

    # Only one call is allowed for ALT
    $$x{ALT}=[$$x{ALT}[0]];
    $$x{FILTER}=[$$x{FILTER}[0]];

    # We are just looking at homozygotes here. Strip
    # out the heterozygote notation.
    $$x{gtypes}{$samplename}{GT}=substr($$x{gtypes}{$samplename}{GT},0,1);

    # It's a SNP if varscan filled in ALT
    my $is_snp=0;
    if($$x{ALT}[0] ne '.'){
      $is_snp=1;
    }

    # some info fields belong in format fields (deleted in the next step from INFO)
    for(qw(ADP)){
      $$x{gtypes}{$samplename}{$_}=$$x{INFO}{$_};
    }
    # I just don't think some fields even belong in the info field
    for(qw(ADP HET NC HOM WT)){
      delete($$x{INFO}{$_});
    }

    # Some FORMAT fields need to be added
    push(@{$$x{FORMAT}}, qw(ADP AF RF FT));

    # Put in allele frequency for the one accepted ALT
    # Alternate calls are in the AD field for varscan
    $$x{gtypes}{$samplename}{AD}||=0;
    $$x{gtypes}{$samplename}{RD}||=0;
    $$x{gtypes}{$samplename}{AF}||=0;
    $$x{gtypes}{$samplename}{RF}||=0;
    if($$x{gtypes}{$samplename}{ADP} > 0){
      $$x{gtypes}{$samplename}{AF}=$$x{gtypes}{$samplename}{AD}/$$x{gtypes}{$samplename}{ADP};
      $$x{gtypes}{$samplename}{RF}=$$x{gtypes}{$samplename}{RD}/$$x{gtypes}{$samplename}{ADP};
    }
    # round the frequencies
    $_=sprintf("%0.2f",$_) for($$x{gtypes}{$samplename}{AF},$$x{gtypes}{$samplename}{RF});

    ###########
    # MASKING
    # Mask low coverage
    $$x{gtypes}{$samplename}{DP} ||= 0;
    if( $$x{gtypes}{$samplename}{DP} < $$settings{coverage} ){
      vcf_markAmbiguous($x,$samplename,$fail_lowCoverage,$settings);
    }

    # Mask indels
    if(length($$x{ALT}[0]) > 1 || length($$x{REF}) > 1){
      vcf_markAmbiguous($x,$samplename,$fail_indel,$settings);
    }
    # Mask sites that are outright masked
    if($$excludeSites{$$x{CHROM}}{$pos}){
      vcf_markAmbiguous($x,$samplename,$fail_masked,$settings);
    }

    # Mask low frequency
    if($is_snp && $$x{gtypes}{$samplename}{AF} < $$settings{altFreq}){
      vcf_markAmbiguous($x,$samplename,$fail_lowRefFreq,$settings);
    }
    if(!$is_snp && $$x{gtypes}{$samplename}{RF} < $$settings{altFreq}){
      vcf_markAmbiguous($x,$samplename,$fail_lowAltFreq,$settings);
    }
    # END MASKING
    ###############
    
    # After all this masking or not masking, 
    # the FT field should match the FILTER field.
    $$x{gtypes}{$samplename}{FT}=join(";",@{$$x{FILTER}});

    print $vcfOut $vcf->format_line($x);
  }
  close $vcfOut;

  return $newVcf;
}

# VCF line manipulations
sub vcf_markAmbiguous{
  my($x,$samplename,$reason,$settings)=@_;

  # Empty the filter if it only says that it passes,
  # because it does not pass anymore.
  @{$$x{FILTER}}=() if($$x{FILTER}[0] eq "PASS");
  # Enter the reason why it fails
  push(@{$$x{FILTER}},$reason);

  # The site failed: GT is unknown. Mark it with
  # a dot.
  $$x{gtypes}{$samplename}{GT}=".";
}

# Get a hash of seqname->{pos} from a bed file
sub readBed{
  my($bed,$settings)=@_;

  my %bed;
  my $bedin=Bio::FeatureIO->new(-format=>"bed",-file=>$bed);
  while(my $feat=$bedin->next_feature){
    my $seqname=$feat->seq_id;
    my $start=$feat->start;
    my $end=$feat->end;
    for my $pos($start..$end){
      $bed{$seqname}{$pos}=1;
    }
  }
  $bedin->close;
  return \%bed;
}

sub usage{
  local $0=fileparse $0;
  "$0: find SNPs using varscan
  Usage: $0 file.bam --reference ref.fasta > file.vcf
  --reference          The reference fasta file
  --numcpus  1         How many cpus to use
  --tempdir  tmp/      A temporary directory to store files
  --coverage 10        Min coverage
  --altFreq  0.75      Min consensus agreement for a SNP
  --region   file.bed  File of positions to include (default with no bed file: read all positions)
  --exclude  file.bed  File of positions to mask    (conflicts with --region; default: don't exclude)
  "
}
