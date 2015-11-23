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

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/@fastaExt @fastqExt @bamExt/;
use Vcf;

my $samplename="SampleNameUnknown"; # to be changed later
$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
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
  $samplename=basename($bam,@bamExt);
  $samplename=~s/\-$refname$//;
  $samplename=basename($samplename,@fastqExt,@fastaExt);
  $samplename=basename($samplename,@fastaExt);

  my $vcf=varscan2($bam,$reference,$settings);
  #backfillVcfValues($vcf,$settings);
  system("cat $vcf");
  die "ERROR: missing tmp vcf file $vcf, or I could not read it" if $?;

  system("rm -f $vcf");

  return 0;
}

sub varscan2{
  my($bam,$reference,$settings)=@_;
  my $b=fileparse($bam);
  my $pileup="$$settings{tempdir}/$b.mpileup";
  return $pileup if(-e $pileup && -s $pileup > 0);
  logmsg "Creating a pileup $pileup";

  # Make sure this whole operation is self-contained
  my $tmpdir="$$settings{tempdir}/$samplename";
  system("rm -rf $tmpdir");
  mkdir($tmpdir);

  # Figure out any regions
  my $regions=`makeRegions.pl --numcpus $$settings{numcpus} --numchunks $$settings{numcpus} $bam`;

  # Figure out mpileup options
  my $xopts="";
  $xopts.="-Q 0 -B "; # assume that all reads have been properly filtered at this point and that mappings are good
  $xopts.="--positions $$settings{region} " if($$settings{region});

  # Make sure that no other pileup file gets in the way.
  system("rm -f $$settings{tempdir}/$b*.mpileup $$settings{tempdir}/$b*.vcf.gz* >&2");

  # Multithread a pileup so that each region gets piped into a file.
  # Then, these individual files can be combined at a later step.
  my $command=qq(echo "$regions" | xargs -P $$settings{numcpus} -n 1 -I {} bash -c '
    echo "MPileup on {}"; 
    pileup=$tmpdir/\$\$.mpileup;
    vcf=$tmpdir/\$\$.vcf;

    # Avoid disk I/O problems.
    sleep \$((\$RANDOM % $$settings{numcpus}))

    samtools mpileup -f $reference $xopts --region "{}" $bam > \$pileup;
    if [ \$? -gt 0 ]; then exit 1; fi;
    size=\$(wc -c < \$pileup);
    if [ "\$size" -lt 1 ]; then
      echo \$pileup is empty. I will not run varscan on it.
      rm -v \$pileup;
      exit 0;
    else
      echo \$pileup is greater than 0 bytes. I will run varscan on it.
    fi;
    varscan.sh mpileup2cns \$pileup --min-coverage $$settings{coverage} --min-var-freq $$settings{altFreq} --output-vcf 1 --min-avg-qual 0 | \\
    sed "s/Sample1/$samplename/" | \\
    bgzip -c > \$vcf.gz.tmp && \
    mv -v \$vcf.gz.tmp \$vcf.gz && \
    tabix \$vcf.gz
  ' >&2);
  logmsg "Running mpileup and varscan:\n  $command";
  system($command);
  die "ERROR with xargs and samtools mpileup" if $?;

  # bcftools concat $TEMPDIR/merged.*.vcf.gz | vcf-sort > $TEMPDIR/concat.vcf
  system("bcftools concat $tmpdir/*.vcf.gz  > $$settings{tempdir}/$b.merged.vcf");
  die if $?;

  # Make sure everything is cleaned up whenever the script ends
  system("rm -rf $tmpdir");

  #backfillVcfValues("$$settings{tempdir}/merged.vcf",$settings);

  return "$$settings{tempdir}/$b.merged.vcf";
}

sub backfillVcfValues{
  my($vcfFile,$bam,$settings)=@_;

  my $excludeSites=readBed($$settings{exclude},$settings) if($$settings{exclude});

  logmsg "Backfilling values in $vcfFile and printing to stdout";
  my $vcf=Vcf->new(file=>$vcfFile);

  # Reasons a SNP doesn't pass
  my $fail_lowCoverage="DP$$settings{coverage}";
  my $fail_lowRefFreq ="RF$$settings{altFreq}";
  my $fail_lowAltFreq ="AF$$settings{altFreq}";
  my $fail_indel      ="isIndel";
  my $fail_masked     ="masked";

  # How was the bam generated?
  my $mapper="unknown";
  my $CL="unknown";
  open(SAMTOOLS,"samtools view -H '$bam' | ") or die "ERROR: could not open $bam with samtools: $!";
  while(<SAMTOOLS>){
    next if(!/^\@PG/);
    chomp;
    my @F=split /\t/;
    for(@F){
      my($key,$value)=split /:/;
      $mapper=$value if($key eq "ID");
      $CL=$value if($key eq "CL");
    }
  }
  
  # Add new headers
  $vcf->add_header_line({key=>'reference',value=>$$settings{reference}});
  $vcf->add_header_line({key=>'mapper',value=>$mapper});
  $vcf->add_header_line({key=>'mapperCL',value=>$CL});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowCoverage, Description=>"Depth is less than $$settings{coverage}, the user-set coverage threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowRefFreq, Description=>"Reference variant consensus is less than $$settings{altFreq}, the user-set threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_lowAltFreq, Description=>"Allele variant consensus is less than $$settings{altFreq}, the user-set threshold"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_indel, Description=>"Indels are not used for analysis in Lyve-SET"});
  $vcf->add_header_line({key=>'FILTER', ID=>$fail_masked, Description=>"This site was masked using a bed file or other means"});

  # Genotype fields
  $vcf->add_header_line({key=>'FORMAT', ID=>'AF', Number=>'A', Type=>'Float', Description=>"Allele Frequency"});
  $vcf->add_header_line({key=>'FORMAT', ID=>'RF', Number=>'A', Type=>'Float', Description=>"Reference Frequency"});
  $vcf->add_header_line({key=>'FORMAT', ID=>'ADP', Number=>'A', Type=>'Integer', Description=>"Depth of bases with Phred score >= 15"});

  # Remove unwanted headers
  $vcf->remove_header_line(key=>'INFO',ID=>'ADP');

  # Figure out the coverage depth
  my $depth=covDepth($bam,$settings);

  # Done with headers; parse them
  $vcf->parse_header();

  # start printing
  print $vcf->format_header();
  while(my $x=$vcf->next_data_hash()){
    my $posId=$$x{CHROM}.':'.$$x{POS};
    my $pos=$$x{POS};
    $$x{gtypes}{$samplename}{DP}=$$depth{$posId} || 0;

    # Only one call is allowed for ALT
    $$x{ALT}=[$$x{ALT}[0]];
    $$x{FILTER}=[$$x{FILTER}[0]];

    # It's a SNP if it isn't the same as REF and isn't a dot
    my $is_snp=0;
    if($$x{ALT}[0] ne $$x{REF} && $$x{ALT}[0] ne '.'){
      $is_snp=1;
    }

    # some info fields belong in format fields (deleted in the next step from INFO)
    for(qw(ADP)){
      $$x{gtypes}{$samplename}{$_}=$$x{INFO}{$_};
    }
    # I just don't think some fields even belong here
    for(qw(ADP HET NC HOM WT)){
      delete($$x{INFO}{$_});
    }

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
    $_=sprintf("%0.3f",$_) for($$x{gtypes}{$samplename}{AF},$$x{gtypes}{$samplename}{RF});

    ###########
    # MASKING
    # Mask low coverage
    if( $$depth{$posId} < $$settings{coverage} ){
      vcf_markAmbiguous($x,$fail_lowCoverage,$settings);
    }

    # Mask indels
    if(length($$x{ALT}[0]) > 1 || length($$x{REF}) > 1){
      vcf_markAmbiguous($x,$fail_indel,$settings);
    }
    # Mask sites that are outright masked
    if($$excludeSites{$$x{CHROM}}{$pos}){
      vcf_markAmbiguous($x,$fail_masked,$settings);
    }

    # Mask low frequency
    if($is_snp && $$x{gtypes}{$samplename}{AF} < $$settings{altFreq}){
      vcf_markAmbiguous($x,$fail_lowRefFreq,$settings);
    }
    if(!$is_snp && $$x{gtypes}{$samplename}{RF} < $$settings{altFreq}){
      vcf_markAmbiguous($x,$fail_lowAltFreq,$settings);
    }
    # END MASKING
    ###############

    print $vcf->format_line($x);
  }

  return 1;
}

sub covDepth{
  my($bam,$settings)=@_;
  my $depthFile="$bam.depth";
  my %depth;
  die if $?;
  open(IN,"samtools depth '$bam' | ") or die "Could not open $bam in samtools:$!";
  while(<IN>){
    my($rseq,$pos,$depth)=split(/\t/);
    chomp($depth);
    $depth{$rseq.":".$pos}=$depth;
  }
  close IN;
  return \%depth;
}

# VCF line manipulations
sub vcf_markAmbiguous{
  my($x,$reason,$settings)=@_;
  my $ucRef=uc($$x{REF}); # uppercase reference, to make it easier for str comparison

  # Figure out the alt allele (N) and the correct genotype integer
  my $gtInt=0; # for the GT, e.g., 1/1 or 2/2
  # Mark if it is the reference allele
  if($ucRef eq $$x{ALT}[0] || $ucRef eq '.'){
    #shift(@{ $$x{ALT} }); # unnecessary to mark the ref base in ALT
    $gtInt=0;
  } else {
    $gtInt=1;
  }
  #push(@{ $$x{ALT} }, "N");
  $$x{ALT}=["N"];
  $$x{gtypes}{$samplename}{GT}="$gtInt/$gtInt";

  # empty the filter and then put on the fail reason
  @{$$x{FILTER}}=() if($$x{FILTER}[0] eq "PASS");
  push(@{$$x{FILTER}},$reason);
  #$$x{FILTER}=[$reason];
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
  --numcpus  1         How many cpus to use
  --tempdir  tmp/      A temporary directory to store files
  --coverage 10        Min coverage
  --altFreq  0.75      Min consensus agreement for a SNP
  --region   file.bed  File of positions to include (default with no bed file: read all positions)
  --exclude  file.bed  File of positions to mask    (conflicts with --region; default: don't exclude)
  "
}
