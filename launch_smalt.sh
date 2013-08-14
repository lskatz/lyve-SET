#!/bin/sh

# qsub -N 'q$b' -cwd -V -o $log/$b.out -e $log/$b.out $scriptsdir/launch_smalt.sh $ref $fastq $bamdir/$b.bam $tmpdir

ref=$1
query=$2
out=$3
tempdir=$4

if [ "$tempdir"="" ]; then
  tempdir="tmp"
  echo "No tempdir was given. Setting it as $tempdir";
fi;

numCpus=24

export ref; # for bioperl, below
b=`basename $out .bam`;
d=`dirname $out`;

if [ "$out" = "" ]; then
  echo usage: $0 ref.fasta query.fastq out.bam tempdir/
  echo "  " The final file will be out.sorted.bam and out.sorted.bam.depth
  exit 1;
fi;

echo "Trimming, removing duplicate reads, and cleaning"
b=`basename $query`;
prefix="$tempdir/$b";
if [ ! -e "$prefix.cleaned.fastq" ]; then
  run_assembly_trimLowQualEnds.pl -c 99999 $query > "$prefix.trimmed.fastq"
  if [ "$?" -gt 0 ]; then exit 1; fi;
  run_assembly_removeDuplicateReads.pl "$prefix.trimmed.fastq" > "$prefix.nodupes.fastq";
  if [ "$?" -gt 0 ]; then exit 1; fi;
  run_assembly_trimClean.pl -i "$prefix.nodupes.fastq" -o "$prefix.cleaned.fastq" --min_length 36
  if [ "$?" -gt 0 ]; then exit 1; fi;
  rm -v "$prefix.nodupes.fastq" "$prefix.trimmed.fastq"
fi;

echo
echo "Mapping against the reference"
if [ ! -e "$b.sorted.bam.depth" ]; then
  if [ `grep '>' $ref | grep -c _` -gt 0 ]; then
    echo "ERROR: You cannot have deflines with underscores in them because of the way this script handles deflines internally."
    echo "Suggestion: change all underscores to dashes"
    echo "  sed -ip 's/_/-/g' $ref # an example "
    exit 1;
  fi
  tmpOut="$out.$RANDOM.tmp";
  echo "$out not found. I'm really doing the mapping now!"
  echo "  Outfile: $tmpOut"
  smalt map -r -1 -f samsoft -n $numCpus $ref "$prefix.cleaned.fastq" | samtools view -bS -T $ref - > $tmpOut;
  if [ "$?" -gt 0 ]; then exit 1; fi;

  echo "Transforming the output file with samtools"
  # sort, index, depth
  sPrefix="$tmpOut.sorted";
  sorted="$sPrefix.bam" # samtools adds the bam later, so I want to keep track of it

  echo "Sorting the temporary bam file into $sorted"
  samtools sort $tmpOut $sPrefix
  samtools index $sorted

  # get the depth and include zero depth
  # http://www.biostars.org/p/19175/
  samtools depth $sorted |\
    perl -MBio::Perl -e 'while(<>){chomp; @F=split /\t/;$depth{$F[0]}{$F[1]}=$F[2];} $in=Bio::SeqIO->new(-file=>"'$ref'");while($seq=$in->next_seq){$max{$seq->id}=$seq->length;} for $chr(keys(%max)){for $i(1..$max{$chr}){$depth{chr}{$i}||=0; print join("\t",$chr,$i,$depth{$chr}{$i})."\n";}}'\
    > "$sorted.depth"

fi

# cleanup
mv -v $sorted $out
rm -v $tmpOut

