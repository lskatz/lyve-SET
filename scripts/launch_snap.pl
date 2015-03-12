#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/fileparse dirname basename/;
use Bio::Perl;
use File::Copy qw/copy move/;

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
  $$settings{samtoolsxopts}||="";

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

  # PE reads or not?  Mapping differently for different types.
  if(is_fastqPE($query)){
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

#    # SNAP checks the read IDs and is incorrectly reading them as wrongly paired.
#    # Therefore I need to make them more kosher with new IDs.
#    # Unfortunately this removes the original information but I am not sure how to
#    # retain it properly while removing information that trips up snap.
#    my $filecount=0;
#    for my $file("$prefix.1.fastq","$prefix.2.fastq"){
#      $filecount++;  # this will be used as the /1 or /2 in the identifier
#      move($file,"$file.tmp");
#      open(FASTQMODDED,">",$file) or die "ERROR: could not open for writing: $file: $!";
#      open(FASTQ,"$file.tmp") or die "ERROR: could not open for reading $file.tmp: $!";
#      my $i=0;
#      my $identifier;
#      while(<FASTQ>){
#        my $mod=++$i % 4;
#        if($mod==1){
#          $identifier="\@$i/$filecount";
#          s/.*/$identifier/;  #just replace the defline with whatever as long as it's consistent
#        } elsif($mod==3){
#          s/.*/\+/;
#        }
#        print FASTQMODDED $_;
#      }
#      close FASTQ;
#      close FASTQMODDED;
#    }

    # mapping ###
    # SNAP gives weird permissions to the output file. Avoid it by
    # creating the file yourself first.
    my $snap=`which snap`; chomp($snap);
    system("touch $tmpOut");
    die if $?;
    my $command="$snap paired $ref.snap '$prefix.1.fastq.gz' '$prefix.2.fastq.gz' -t $$settings{numcpus} -so -o $tmpOut -h 1 -C++ --b -I";
    system($command);
    die "ERROR with command: $!\n $command" if $?;

    system("rm -v '$prefix.1.fastq.gz' '$prefix.2.fastq.gz'"); die if $?;
  } else {
    ##########
    # SE reads
    ##########
    system("gunzip -c '$query' > $prefix.SE.fastq");
    die "Problem with gunzip" if $?;
    system("touch $tmpOut $tmpOut.bai");
    die if $?;
    system("snap single $ref.snap '$prefix.SE.fastq' -t $$settings{numcpus} -so -o $tmpOut -h 1");
    if($?){
      # longer reads will cause an error.  Try the snapxl binary instead just in case
      logmsg "Snap failed with an error.  Trying snapxl";
      system("snapxl single $ref.snap '$prefix.SE.fastq' -t $$settings{numcpus} -so -o $tmpOut");
      die "ERROR: snpxl failed. Maybe you don't have it in your path? It can be compiled from snap using 'make snapxl' from the snap package (https://github.com/amplab/snap)" if $?;
    }

    system("rm -v '$prefix.SE.fastq'"); die if $?;
  }

  # SNAP includes the entire defline for some reason, but the 
  # part in the defline after the space needs to be removed for 
  # compatibility with samtools
  my $seqin=Bio::SeqIO->new(-file=>$ref);
  while(my $seq=$seqin->next_seq){
    my $desc=$seq->desc;
    $desc=~s/\s/_/g;
    logmsg "Removing _$desc from sam file";
    system("sed -i 's/_$desc//' $tmpOut");
    logmsg "WARNING: could not replace a desc from the defline: $!\n  The temporary sam file might not be compatible with the rest of the automation." if $?;
  }


  # Do some post-mapping filtering
  $$settings{samtoolsxopts}.="-F 4 ";
  my $regions="$$settings{refdir}/unmaskedRegions.bed";
  if(-e $regions && -s $regions > 0){
    logmsg "Found bed file of regions to accept and so I am using it! $regions";
    $$settings{samtoolsxopts}.="-L $regions ";
  }
  system("mv -v $tmpOut $tmpOut.unfiltered && samtools view $$settings{samtoolsxopts} -bSh -T $ref $tmpOut.unfiltered > $tmpOut.bam");
  die "ERROR with samtools view" if $?;
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
  "
}

