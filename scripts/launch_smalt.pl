#!/usr/bin/env perl

# qsub -N 'q$b' -cwd -V -o $log/$b.out -e $log/$b.out $scriptsdir/launch_smalt.sh $ref $fastq $bamdir/$b.bam $tmpdir

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/fileparse dirname basename/;
use File::Copy qw/cp/;
use Bio::Perl;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use LyveSET qw/@fastqExt logmsg/;

$0=fileparse $0;
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help reference=s fastq=s bam=s tempdir=s numcpus=i smaltxopts=s pairedend=i minPercentIdentity=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||="tmp";
  $$settings{pairedend}||=0;
  $$settings{minPercentIdentity}||=95;
  $$settings{minPercentIdentity}=sprintf("%0.2f",$$settings{minPercentIdentity}/100);

  # smalt extra options
  $$settings{smaltxopts}||="";
  #$$settings{smaltxopts}.=" -i 1000 -y 0.95 -f samsoft -n $$settings{numcpus}";
  $$settings{smaltxopts}.=" -i 1000 -y $$settings{minPercentIdentity} -r -1 -f samsoft -n $$settings{numcpus}";
  $$settings{samtoolsxopts}||="";

  for(qw(reference fastq bam)){
    die "ERROR: need option $_\n".usage($settings) if(!$$settings{$_});
  }
  die usage($settings) if($$settings{help});

  mkdir $$settings{tempdir} if(!-d $$settings{tempdir});

  my $fastq=$$settings{fastq};
  my $bam=$$settings{bam};
  my $reference=$$settings{reference};
  $$settings{refdir}||=dirname($reference);

  # Check if the reference genome was indexed.
  # I think it has to die and not index the genome just in case other
  # instances of this script are trying to do the same thing.
  for my $indexFile("$reference.sma","$reference.smi"){
    die "ERROR: Could not find $indexFile, which makes me think that the reference genome has not been indexed with `smalt index`." if(!-e $indexFile);
  }

  mapReads($fastq,$bam,$reference,$settings);

  return 0;
}

sub mapReads{
  my($query,$bam,$ref,$settings)=@_;
  if(-s "$bam"){
    logmsg "Found $bam\n  I will not re-map.";
    return 1;
  }

  my ($b,$infilePath,$infileExt)=fileparse($query,@fastqExt);
  my $prefix="$$settings{tempdir}/$b";
  
  my $RANDOM=rand(999999);
  my $tmpOut="$bam.$RANDOM.tmp";
  my $tmpSamOut="$prefix.tmp.sam";
  logmsg "$bam not found. I'm really doing the mapping now!\n  Outfile: $tmpOut";

  # Find out if there is a set of coordinates to use
  my $regions="$$settings{refdir}/unmaskedRegions.bed";
  if(-e $regions && -s $regions > 0){
    logmsg "Found bed file of regions to accept and so I am using it! $regions";
    $$settings{samtoolsxopts}.="-L $regions ";
  }
  #$$settings{samtoolsxopts}.="-F 4 "; # only keep mapped reads to save on space

  # PE reads or not?  Mapping differently for different types.
  if(is_fastqPE($query,$settings)){
    # deshuffle the reads
    logmsg "Deshuffling to $prefix.1.fastq and $prefix.2.fastq";
    system("run_assembly_shuffleReads.pl -d '$query' 1>'$prefix.1.fastq' 2>'$prefix.2.fastq'");
    die "Problem with deshuffling reads! I am assuming paired end reads." if $?;

    # Mapping
    system("smalt map $$settings{smaltxopts} $ref '$prefix.1.fastq' '$prefix.2.fastq' > $tmpSamOut");
    die if $?;
    system("rm -v '$prefix.1.fastq' '$prefix.2.fastq'"); die if $?;
  } else {
    # remove any paired end parameters that would cause an error
    $$settings{smaltxopts}=~s/-i\s+\d+//;

    # Need to think about gzip vs un-gzipped files
    if($infileExt=~/\.gz$/){
      system("gunzip -c '$query' > $prefix.SE.fastq");
      die "Problem gunzipping $query to $prefix.SE.fastq" if $?;
    } else {
      cp($query,"$prefix.SE.fastq")
        or die "ERROR: could not copy to $prefix.SE.fastq";
    }

    # Mapping
    system("smalt map $$settings{smaltxopts} $ref '$prefix.SE.fastq' > $tmpSamOut"); 
    die if $?;
    system("rm -v '$prefix.SE.fastq'"); die if $?;
  }

  # Convert to bam
  system("samtools view $$settings{samtoolsxopts} -bS -T $ref $tmpSamOut > $tmpOut");
  die "ERROR with samtools view $$settings{samtoolsxopts}" if $?;

  logmsg "Transforming the output file with samtools";
  my $sPrefix="$tmpOut.sorted";
  my $sorted="$sPrefix.bam"; # samtools adds the bam later, so I want to keep track of it

  logmsg "Sorting the temporary bam file into $sorted";
  system("samtools sort $tmpOut $sPrefix 2>&1"); die "ERROR with 'samtools sort $tmpOut $sPrefix'" if $?;
  system("samtools index $sorted 2>&1"); die "ERROR with 'samtools index $sorted'" if $?;

  system("set_samtools_depth.pl $sorted 2>&1");
  die if $?;

  # cleanup
  system("mv -v $sorted $bam"); die if $?;
  system("mv -v $sorted.bai $bam.bai"); die if $?;
  unlink("$sorted.depth.gz");
  system("rm -v $tmpOut $tmpSamOut"); die if $?;

  return 1;
}

# Use CGP to determine if a file is PE or not
# settings:  checkFirst is an integer to check the first X deflines
sub is_fastqPE($;$){
  my($fastq,$settings)=@_;

  # If PE was specified, return true for the value 2 (PE)
  # and 0 for the value 1 (SE)
  if($$settings{pairedend}){
    return ($$settings{pairedend}-1);
  }

  # if checkFirst is undef or 0, this will cause it to check at least the first 50 entries.
  # 50 reads is probably enough to make sure that it's shuffled (1/2^25 chance I'm wrong)
  $$settings{checkFirst}||=50;
  $$settings{checkFirst}=50 if($$settings{checkFirst}<2);
  my $is_PE=`run_assembly_isFastqPE.pl '$fastq' --checkfirst $$settings{checkFirst}`;
  chomp($is_PE);
  return $is_PE;
}


sub usage{
  my($settings)=@_;
  "Maps a read set against a reference genome using smalt. Output file will be file.bam and file.bam.depth
  Usage: $0 -f file.fastq -b file.bam -t tmp/ -r reference.fasta
  -t tmp to set the temporary directory as 'tmp'
  --numcpus 1 number of cpus to use
  -s '' Extra smalt map options (not validated). Default: $$settings{smaltxopts} 
  --pairedend <0|1|2> For 'auto', single-end, or paired-end respectively. Default: auto (0).
  --minPercentIdentity 95  The percent identity between a read and its match before it can be mapped
  "
}

