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
  GetOptions($settings,qw(help reference=s fastq=s bam=s tempdir=s clean! numcpus=i));
  for(qw(reference fastq bam)){
    die "ERROR: need option $_\n".usage() if(!$$settings{$_});
  }
  $$settings{numcpus}||=1;
  $$settings{tempdir}||="tmp";

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
    system("smalt map -r -1 -f samsoft -n $$settings{numcpus} $ref '$prefix.1.fastq' '$prefix.2.fastq' | samtools view -bS -T $ref - > $tmpOut");
    die if $?;
    system("rm -v '$prefix.1.fastq' '$prefix.2.fastq'"); die if $?;
  } else {
    system("gunzip -c '$query' > $prefix.SE.fastq");
    system("smalt map -r -1 -f samsoft -n $$settings{numcpus} $ref '$prefix.SE.fastq' | samtools view -bS -T $ref - > $tmpOut");
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

# taken from AKUtils
# See whether a fastq file is paired end or not. It must be in a velvet-style shuffled format.
# In other words, the left and right sides of a pair follow each other in the file.
# params: fastq file and settings
# fastq file can be gzip'd
# settings:  checkFirst is an integer to check the first X deflines
sub is_fastqPE($;$){
  my($fastq,$settings)=@_;

  # if checkFirst is undef or 0, this will cause it to check at least the first 20 entries.
  # 20 reads is probably enough to make sure that it's shuffled (1/2^10 chance I'm wrong)
  $$settings{checkFirst}||=20;
  $$settings{checkFirst}=20 if($$settings{checkFirst}<2);

  # get the deflines
  my @defline;
  my $numEntries=0;
  my $i=0;
  my $fp;
  if($fastq=~/\.gz$/){
    open($fp,"gunzip -c '$fastq' |") or die "Could not open $fastq for reading: $!";
  }else{
    open($fp,"<",$fastq) or die "Could not open $fastq for reading: $!";
  }
  my $discard;
  while(my $defline=<$fp>){
    next if($i++ % 4 != 0);
    chomp($defline);
    $defline=~s/^@//;
    push(@defline,$defline);
    $numEntries++;
    last if($numEntries > $$settings{checkFirst});
  }
  close $fp;

  # it is paired end if it validates with any naming system
  my $is_pairedEnd=_is_fastqPESra(\@defline,$settings) || _is_fastqPECasava18(\@defline,$settings) || _is_fastqPECasava17(\@defline,$settings);

  return $is_pairedEnd;
}

sub _is_fastqPESra{
  my($defline,$settings)=@_;
  my @defline=@$defline; # don't overwrite $defline by mistake
  
  for(my $i=0;$i<@defline-1;$i++){
    my($genome,$info1,$info2)=split(/\s+/,$defline[$i]);
    if(!$info2){
      return 0;
    }
    my($instrument,$flowcellid,$lane,$x,$y,$X,$Y)=split(/:/,$info1);
    my($genome2,$info3,$info4)=split(/\s+/,$defline[$i+1]);
    my($instrument2,$flowcellid2,$lane2,$x2,$y2,$X2,$Y2)=split(/:/,$info3);
    $_||="" for($X,$Y,$X2,$Y2); # these variables might not be present
    if($instrument ne $instrument2 || $flowcellid ne $flowcellid2 || $lane ne $lane2 || $x ne $x2 || $y ne $y2 || $X ne $X2 || $Y ne $Y2){
      return 0;
    }
  }
  return 1;
}

sub _is_fastqPECasava18{
  my($defline,$settings)=@_;
  my @defline=@$defline;

  for(my $i=0;$i<@defline-1;$i++){
    my($instrument,$runid,$flowcellid,$lane,$tile,$x,$yandmember,$is_failedRead,$controlBits,$indexSequence)=split(/:/,$defline[$i]);
    my($y,$member)=split(/\s+/,$yandmember);

    my($inst2,$runid2,$fcid2,$lane2,$tile2,$x2,$yandmember2,$is_failedRead2,$controlBits2,$indexSequence2)=split(/:/,$defline[$i+1]);
    my($y2,$member2)=split(/\s+/,$yandmember2);

    # Instrument, etc must be the same.
    # The member should be different, usually "1" and "2"
    if($instrument ne $inst2 || $runid ne $runid2 || $flowcellid ne $fcid2 || $tile ne $tile2 || $member>=$member2){
      return 0;
    }
  }
  return 1;
}

# This format is basically whether the ends of the defline alternate 1 and 2.
sub _is_fastqPECasava17{
  my($defline,$settings)=@_;
  my @defline=@$defline;
  for(my $i=0;$i<@defline-1;$i++){
    # Get each member number but return false if it doesn't even exist.
    my ($member1,$member2);
    if($defline[$i] =~ m/(\d+)$/){
      $member1=$1;
    } else {
      return 0;
    }
    if($defline[$i+1] =~ /(\d+)$/){
      $member2=$1;
    } else {
      return 0;
    }

    # The test is whether member1 is less than member2.
    # They can't be equal either.
    if($member1 >= $member2){
      return 0;
    }
  }

  return 1;
}


sub usage{
  "Maps a read set against a reference genome using smalt. Output file will be file.bam and file.bam.depth
  Usage: $0 -f file.fastq -b file.bam -t tmp/ -r reference.fasta
  --noclean if you don't want to clean the reads internally
  -t tmp to set the temporary directory as 'tmp'
  --numcpus 1 number of cpus to use
  "
}

