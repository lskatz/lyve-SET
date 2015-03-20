#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/fileparse dirname basename/;
use Bio::Perl;
use File::Copy qw/copy move/;

$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={clean=>0};
  GetOptions($settings,qw(help reference=s fastq=s bam=s tempdir=s clean! numcpus=i pairedend=i));
  for(qw(reference fastq bam)){
    die "ERROR: need option $_\n".usage() if(!$$settings{$_});
  }
  $$settings{numcpus}||=1;
  $$settings{tempdir}||="tmp";
  $$settings{samtoolsxopts}||="";
  $$settings{pairedend}||=0;

  mkdir $$settings{tempdir} if(!-d $$settings{tempdir});
  logmsg "Temporary directory is $$settings{tempdir}";

  my $fastq=$$settings{fastq};
  my $bam=$$settings{bam};
  my $reference=$$settings{reference};
  $$settings{refdir}||=dirname($reference);

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
  
  my $tmpOut="$bam.sam";
  logmsg "$bam not found. I'm really doing the mapping now!\n  Outfile: $tmpOut";

  my $snap=`which snap`; chomp($snap);
  # PE reads or not?  Mapping differently for different types.
  if(is_fastqPE($query,$settings)){
    ###########
    # PE reads
    # #########

    # deshuffle the reads
    logmsg "Deshuffling to $prefix.1.fastq and $prefix.2.fastq";
    system("run_assembly_shuffleReads.pl -d '$query' 1>'$prefix.1.fastq' 2>'$prefix.2.fastq'");
    die "Problem with deshuffling reads! I am assuming paired end reads." if $?;
    # SNAP 
    system("gzip -vf1 $prefix.1.fastq $prefix.2.fastq");
    die "ERROR with gzip" if $?;

    # mapping ###
    # SNAP gives weird permissions to the output file. Avoid it by
    # creating the file yourself first.
    system("touch $tmpOut");
    die if $?;
    my $snap_command="$snap paired $ref.snap '$prefix.1.fastq.gz' '$prefix.2.fastq.gz' -t $$settings{numcpus} -so -o $tmpOut -x -f -C++ --b -I";
    my $snapxl_command=$snap_command;
    $snapxl_command=~s/snap/snapxl/;
    system("$snap_command || $snapxl_command");
    die "ERROR: snap failed! $!   $snap_command || $snapxl_command" if($?);

    system("rm -v '$prefix.1.fastq.gz' '$prefix.2.fastq.gz'"); die if $?;
  } else {
    ##########
    # SE reads
    ##########
    system("gunzip -c '$query' > $prefix.SE.fastq");
    die "Problem with gunzip" if $?;
    system("touch $tmpOut $tmpOut.bai");
    die if $?;
    my $snap_command="$snap single $ref.snap '$prefix.SE.fastq' -t $$settings{numcpus} -so -o $tmpOut -h 1 --b -I"; # -C parameter was causing segfault
    my $snapxl_command=$snap_command;
    $snapxl_command=~s/snap/snapxl/;

    system("$snap_command || $snapxl_command");
    # try one more time
    if($?){
      system("$snap_command || $snapxl_command");
    }
    die "ERROR: snap failed! $!   $snap_command || $snapxl_command" if($?);

    system("rm -v '$prefix.SE.fastq'"); die if $?;
  }

  ## SNAP includes the entire defline for some reason, but the 
  ## part in the defline after the space needs to be removed for 
  ## compatibility with samtools
  #my $seqin=Bio::SeqIO->new(-file=>$ref);
  #while(my $seq=$seqin->next_seq){
  #  my $desc=$seq->desc;
  #  $desc=~s/\s/_/g;
  #  logmsg "Removing _$desc from sam file";
  #  system("sed -i 's/_$desc//' $tmpOut");
  #  logmsg "WARNING: could not replace a desc from the defline: $!\n  The temporary sam file might not be compatible with the rest of the automation." if $?;
  #}


  # Do some post-mapping filtering
  $$settings{samtoolsxopts}.="-F 4 ";
  my $regions="$$settings{refdir}/unmaskedRegions.bed";
  if(-e $regions && -s $regions > 0){
    logmsg "Found bed file of regions to accept and so I am using it! $regions";
    $$settings{samtoolsxopts}.="-L $regions ";
  }
  my $samtoolsView="mv -v $tmpOut $tmpOut.unfiltered && samtools view $$settings{samtoolsxopts} -bSh -T $ref $tmpOut.unfiltered > $tmpOut.bam";
  system($samtoolsView);
  die "ERROR with samtools view\n  $samtoolsView" if $?;
  unlink("$tmpOut.unfiltered");

  # Sort/index
  my $sorted="$bam.sorted.bam";
  system("samtools sort $tmpOut.bam $bam.sorted"); die "ERROR with samtools sort" if $?;
  system("samtools index $sorted"); die "ERROR with samtools index" if $?;
  system("rm -v $tmpOut.bam"); die "ERROR could not delete $tmpOut" if $?;

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
  system("gzip -9v $bam.depth"); die if $?; # save some space
  #system("rm -v $tmpOut"); die if $?;

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
  "Maps a read set against a reference genome using snap. Output file will be file.bam and file.bam.depth
  Usage: $0 -f file.fastq -b file.bam -t tmp/ -r reference.fasta
  -t tmp to set the temporary directory as 'tmp'
  --numcpus 1 number of cpus to use
  --pairedend <0|1|2> For 'auto', single-end, or paired-end respectively. Default: auto (0).
  "
}

