readDir=$1
bamDir=$2
vcfDir=$3
logDir=$4
msaDir=$5

if [ "$msaDir" = "" ]; then
  echo usage: $0 readDir/ bamDir/ vcfDir/ logDir/ msaDir/
  exit 1;
fi

echo "This is not correctly programmed yet."
exit 1

# map reads
(for i in $readDir/*.fastq; do b=`basename $i .fastq`; qsub -N "q$b" -cwd -V -o log/$b.out -e log/$b.out scripts/launch_smalt.sh reference/2010EL-1786.gbk.fasta $i bam/$b.bam; done;)

# variant calls
(for i in $bamDir/*.sorted.bam; do b=`basename $i .sorted.bam`; qsub -N "q$b" -cwd -V -o $logDir/$b.out -e $logDir/$b.out scripts/launch_freebayes.sh reference/2010EL-1786.gbk.fasta $i $vcfDir/$b.vcf; done;)

# variants to MSA
qsub -cwd -o out -e out -V scripts/launch_vcfToAlignment.sh $bamDir/ $vcfDir/ reference/2010EL-1786.gbk.fasta <(echo "") $msaDir/out.aln.fas

# tree
(cd $msaDir; qsub -cwd -o out -e out -V ../scripts/launch_raxml.sh out.aln.fas out; cd ..)
