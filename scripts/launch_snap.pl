#!/usr/bin/env perl

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
  my $settings={clean=>0};
  GetOptions($settings,qw(help reference=s fastq=s bam=s tempdir=s clean! numcpus=i));
  for(qw(reference fastq bam)){
    die "ERROR: need option $_\n".usage() if(!$$settings{$_});
  }
  $$settings{numcpus}||=1;
  $$settings{tempdir}||="tmp";

  mkdir $$settings{tempdir} if(!-d $$settings{tempdir});
  logmsg "Temporary directory is $$settings{tempdir}";

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
  die if $?;
  system("run_assembly_removeDuplicateReads.pl '$query' > '$prefix.nodupes.fastq'");
  die if $?;
  system("run_assembly_trimClean.pl --numcpus $$settings{numcpus} -i '$prefix.nodupes.fastq' -o '$prefix.cleaned.fastq' --auto");
  die if $?;
  system("rm -v '$prefix.nodupes.fastq' '$prefix.trimmed.fastq'");

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

  my $tmpOut="$bam.tmp.bam";
  logmsg "$bam not found. I'm really doing the mapping now!\n  Outfile: $tmpOut";

  # PE reads or not?  Mapping differently for different types.
  if(is_fastqPE($query)){
    ###########
    # PE reads
    # #########

    # deshuffle the reads
    logmsg "Deshuffling to $prefix.1.fastq and $prefix.2.fastq";
    system("run_assembly_shuffleReads.pl -d '$query' 1>'$prefix.1.fastq' 2>'$prefix.2.fastq'");
    die "Problem with deshuffling reads! I am assuming paired end reads." if $?;

    # mapping ###
    # SNAP gives weird permissions to the output file. Avoid it by
    # creating the file yourself first.
    system("touch $tmpOut $tmpOut.bai");
    die if $?;
    system("snap paired $ref.snap '$prefix.1.fastq' '$prefix.2.fastq' -t $$settings{numcpus} -so -o $tmpOut");
    die if $?;

    system("rm -v '$prefix.1.fastq' '$prefix.2.fastq'"); die if $?;
  } else {
    ##########
    # SE reads
    ##########
    system("gunzip -c '$query' > $prefix.SE.fastq");
    die "Problem with gunzip" if $?;
    system("touch $tmpOut $tmpOut.bai");
    die if $?;
    system("snap single $ref.snap '$prefix.SE.fastq' -t $$settings{numcpus} -so -o $tmpOut");
    if($?){
      # longer reads will cause an error.  Try the snapxl binary instead just in case
      logmsg "Snap failed with an error.  Trying snapxl";
      system("snapxl single $ref.snap '$prefix.SE.fastq' -t $$settings{numcpus} -so -o $tmpOut");
      die "ERROR: snpxl failed. Maybe you don't have it in your path? It can be compiled from snap using 'make snapxl' from the snap package (https://github.com/amplab/snap)" if $?;
    }

    system("rm -v '$prefix.SE.fastq'"); die if $?;
  }
  my $sorted=$tmpOut; # because I am too lazy to rename variables

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
  #system("rm -v $tmpOut"); die if $?;

  return 1;
}

# taken from AKUtils
# See whether a fastq file is paired end or not. It must be in a velvet-style shuffled format.
# In other words, the left and right sides of a pair follow each other in the file.
# params: fastq file and settings
# fastq file can be gzip'd
# settings:  checkFirst is an integer to check the first X deflines
# TODO just extract IDs and send them to the other _sub()
sub is_fastqPE($;$){
  my($fastq,$settings)=@_;

  # if checkFirst is undef or 0, this will cause it to check at least the first 20 entries.
  $$settings{checkFirst}||=20;
  $$settings{checkFirst}=20 if($$settings{checkFirst}<2);

  # it is paired end if it validates with any naming system
  my $is_pairedEnd=_is_fastqPESra($fastq,$settings) || _is_fastqPECasava18($fastq,$settings) || _is_fastqPECasava17($fastq,$settings);

  return $is_pairedEnd;
}

sub _is_fastqPESra{
  my($fastq,$settings)=@_;
  my $numEntriesToCheck=$$settings{checkFirst}||20;

  my $numEntries=0;
  my $fp;
  if($fastq=~/\.gz$/){
    open($fp,"gunzip -c '$fastq' |") or die "Could not open $fastq for reading: $!";
  }else{
    open($fp,"<",$fastq) or die "Could not open $fastq for reading: $!";
  }
  my $discard;
  while(<$fp>){
    chomp;
    s/^@//;
    my($genome,$info1,$info2)=split(/\s+/);
    if(!$info2){
      close $fp;
      return 0;
    }
    my($instrument,$flowcellid,$lane,$x,$y,$X,$Y)=split(/:/,$info1);
    $discard=<$fp> for(1..3); # discard the sequence and quality of the read for these purposes
    my $secondId=<$fp>;
    my($genome2,$info3,$info4)=split(/\s+/,$secondId);
    my($instrument2,$flowcellid2,$lane2,$x2,$y2,$X2,$Y2)=split(/:/,$info3);
    $_||="" for($X,$Y,$X2,$Y2); # these variables might not be present
    if($instrument ne $instrument2 || $flowcellid ne $flowcellid2 || $lane ne $lane2 || $x ne $x2 || $y ne $y2 || $X ne $X2 || $Y ne $Y2){
      close $fp;
      return 0;
    }
    $discard=<$fp> for(1..3);
    $numEntries+=2;
    last if($numEntries > $numEntriesToCheck);
  }
  return 1;
}

sub _is_fastqPECasava18{
  my($fastq,$settings)=@_;
  my $numEntriesToCheck=$$settings{checkFirst}||20;

  my $numEntries=0;
  my $fp;
  if($fastq=~/\.gz$/){
    open($fp,"gunzip -c '$fastq' |") or die "Could not open $fastq for reading: $!";
  }else{
    open($fp,"<",$fastq) or die "Could not open $fastq for reading: $!";
  }
  while(<$fp>){
    chomp;
    s/^@//;
    my($instrument,$runid,$flowcellid,$lane,$tile,$x,$yandmember,$is_failedRead,$controlBits,$indexSequence)=split(/:/,$_);
    my $discard;
    $discard=<$fp> for(1..3); # discard the sequence and quality of the read for these purposes
    my($y,$member)=split(/\s+/,$yandmember);

    # if all information is the same, except the member (1 to 2), then it is still paired until this point.
    my $secondId=<$fp>;
    chomp $secondId;
    $secondId=~s/^@//;
    my($inst2,$runid2,$fcid2,$lane2,$tile2,$x2,$yandmember2,$is_failedRead2,$controlBits2,$indexSequence2)=split(/:/,$secondId);
    $discard=<$fp> for(1..3); # discard the sequence and quality of the read for these purposes
    my($y2,$member2)=split(/\s+/,$yandmember2);

    if($instrument ne $inst2 || $runid ne $runid2 || $flowcellid ne $fcid2 || $tile ne $tile2 || $member!=1 || $member2!=2){
      #logmsg "Failed!\n$instrument,$runid,$flowcellid,$lane,$tile,$x,$yandmember,$is_failedRead,$controlBits,$indexSequence\n$inst2,$runid2,$fcid2,$lane2,$tile2,$x2,$yandmember2,$is_failedRead2,$controlBits2,$indexSequence2\n";
      close $fp;
      return 0;
    }

    $numEntries+=2;
    last if($numEntries>$numEntriesToCheck);
  }

  close $fp;

  return 1;
}

sub _is_fastqPECasava17{
  my($fastq,$settings)=@_;
  # 20 reads is probably enough to make sure that it's shuffled (1/2^20 chance I'm wrong)
  my $numEntriesToCheck=$$settings{checkFirst}||20;
  my $numEntries=0;
  my $fp;
  if($fastq=~/\.gz$/){
    open($fp,"gunzip -c '$fastq' |") or die "Could not open $fastq for reading: $!";
  }else{
    open($fp,"<",$fastq) or die "Could not open $fastq for reading: $!";
  }
  while(my $read1Id=<$fp>){
    my $discard;
    $discard=<$fp> for(1..3);
    my $read2Id=<$fp>;
    $discard=<$fp> for(1..3);

    if($read1Id!~/\/1$/ || $read2Id!~/\/2$/){
      close $fp;
      return 0;
    }

    $numEntries+=2;
    last if($numEntries>=$numEntriesToCheck);
  }
  close $fp;

  return 1;
}


sub usage{
  "Maps a read set against a reference genome using snap. Output file will be file.bam and file.bam.depth
  Usage: $0 -f file.fastq -b file.bam -t tmp/ -r reference.fasta
  -t tmp to set the temporary directory as 'tmp'
  --numcpus 1 number of cpus to use
  "
}

