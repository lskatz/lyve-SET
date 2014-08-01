#!/usr/bin/env perl

# qsub -N 'q$b' -cwd -V -o $log/$b.out -e $log/$b.out $scriptsdir/launch_smalt.sh $ref $fastq $bamdir/$b.bam $tmpdir

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Bio::Perl;

$0=fileparse $0;
sub logmsg{print STDERR "@_\n";}
exit(main());

sub main{
  my $settings={clean=>1};
  GetOptions($settings,qw(help reference=s fastq=s bam=s tempdir=s clean! numcpus=i smaltxopts=s)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{tempdir}||="tmp";

  # smalt extra options
  $$settings{smaltxopts}||="";
  $$settings{smaltxopts}.=" -y 0.95 -r -1 -f samsoft -n $$settings{numcpus}";

  for(qw(reference fastq bam)){
    die "ERROR: need option $_\n".usage($settings) if(!$$settings{$_});
  }
  die usage($settings) if($$settings{help});

  mkdir $$settings{tempdir} if(!-d $$settings{tempdir});

  my $fastq=$$settings{fastq};
  my $bam=$$settings{bam};
  my $reference=$$settings{reference};

  logmsg "Cleaning was disabled. I will not clean the reads before mapping" if(!$$settings{clean});
  $fastq=cleanReads($fastq,$settings) if($$settings{clean});
  mapReads($fastq,$bam,$reference,$settings);

  return 0;
}

sub cleanReads{
  my($query,$settings)=@_;
  logmsg "Trimming, removing duplicate reads, and cleaning";
  my $b=fileparse $query;
  my $prefix="$$settings{tempdir}/$b";
  #system("run_assembly_trimLowQualEnds.pl -c 99999 $query > '$prefix.trimmed.fastq'"); 
  die if $?;
  #system("run_assembly_removeDuplicateReads.pl '$prefix.trimmed.fastq' > '$prefix.nodupes.fastq'");
  system("run_assembly_removeDuplicateReads.pl '$query' > '$prefix.nodupes.fastq'");
  die if $?;
  system("run_assembly_trimClean.pl --numcpus $$settings{numcpus} -i '$prefix.nodupes.fastq' -o '$prefix.cleaned.fastq' --min_length 36");
  die if $?;
  system("rm -v '$prefix.nodupes.fastq' '$prefix.trimmed.fastq'");
  #die if $?;

  return "$prefix.cleaned.fastq";
}

sub mapReads{
  my($query,$bam,$ref,$settings)=@_;
  if(-s "$bam.depth"){
    logmsg "Found $bam.depth\n  I will not re-map.";
    return 1;
  }

  my $b=fileparse $query;
  my $prefix="$$settings{tempdir}/$b";
  
  # don't allow underscores in the reference
  my $numUnderscores=`grep '>' $ref | grep -c _`; chomp($numUnderscores);
  if($numUnderscores > 0){
    logmsg "ERROR: You cannot have deflines with underscores in them because of the way this script handles deflines internally.";
    logmsg "Suggestion: change all underscores to dashes";
    logmsg "  sed -ip 's/_/-/g' $ref # an example ";
    logmsg "TODO: automate this step";
    die;
  }

  my $RANDOM=rand(999999);
  my $tmpOut="$bam.$RANDOM.tmp";
  logmsg "$bam not found. I'm really doing the mapping now!\n  Outfile: $tmpOut";

  # PE reads or not?  Mapping differently for different types.
  if(is_fastqPE($query)){
    # deshuffle the reads
    logmsg "Deshuffling to $prefix.1.fastq and $prefix.2.fastq";
    system("run_assembly_shuffleReads.pl -d '$query' 1>'$prefix.1.fastq' 2>'$prefix.2.fastq'");
    die "Problem with deshuffling reads! I am assuming paired end reads." if $?;

    # mapping
    system("smalt map $$settings{smaltxopts} $ref '$prefix.1.fastq' '$prefix.2.fastq' | samtools view -F 4 -bS -T $ref - > $tmpOut");
    die if $?;
    system("rm -v '$prefix.1.fastq' '$prefix.2.fastq'"); die if $?;
  } else {
    system("gunzip -c '$query' > $prefix.SE.fastq");
    die if $?;
    system("smalt map $$settings{smaltxopts} $ref '$prefix.SE.fastq' | samtools view -bS -T $ref - > $tmpOut");
    die if $?;
    system("rm -v '$prefix.SE.fastq'"); die if $?;
  }

  logmsg "Transforming the output file with samtools";
  my $sPrefix="$tmpOut.sorted";
  my $sorted="$sPrefix.bam"; # samtools adds the bam later, so I want to keep track of it

  logmsg "Sorting the temporary bam file into $sorted";
  system("samtools sort $tmpOut $sPrefix"); die if $?;
  system("samtools index $sorted"); die if $?;

  logmsg "Getting the depth at each position in the sorted bam file";
  # read the depth, but unfortunately it skips over zero-depth sites
  my %depth;
  open(BAMDEPTH,"samtools depth $sorted | ") or die "ERROR: could not open the bam file with samtools: $!";
  while(<BAMDEPTH>){
    chomp;
    my @F=split /\t/;
    $depth{$F[0]}{$F[1]}=$F[2];
  }
  close BAMDEPTH;
 
  # get zero-depth information in there too. First find the total length of each contig.
  open(FIXEDDEPTH,">","$sorted.depth") or die "ERROR: could not write to $sorted.depth: $!";
  my %max;
  my $in=Bio::SeqIO->new(-file=>$ref);
  while(my $seq=$in->next_seq){
    $max{$seq->id}=$seq->length;
  }
  for my $contig(keys(%max)){
    for my $i(1..$max{$contig}){
      $depth{$contig}{$i}||=0; 
      print FIXEDDEPTH join("\t",$contig,$i,$depth{$contig}{$i})."\n";
    }
  }
  close FIXEDDEPTH;
  
  system("mv -v $sorted $bam"); die if $?;
  system("mv -v $sorted.bai $bam.bai"); die if $?;
  system("mv -v $sorted.depth $bam.depth"); die if $?;
  system("rm -v $tmpOut"); die if $?;

  return 1;
}

# Use CGP to determine if a file is PE or not
# settings:  checkFirst is an integer to check the first X deflines
sub is_fastqPE($;$){
  my($fastq,$settings)=@_;

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
  --noclean if you don't want to clean the reads internally
  -t tmp to set the temporary directory as 'tmp'
  --numcpus 1 number of cpus to use
  -s '' Extra smalt map options (not validated). Default: $$settings{smaltxopts} 
  "
}

