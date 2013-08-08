#!/bin/sh

bamDir=$1
vcfDir=$2
ref=$3
out=$4

if [ "$out" = "" ]; then
  echo usage: $0 bam/ vcf/ ref.fasta badSites.txt out.aln.fas
  echo "  bam directory's *.sorted.bam will be read"
  echo "  vcf dir's *.vcf will be read"
  echo "  vcfDir/*.badsites.txt will be used as 'bad sites', as produced by filterVcf.pl"
  exit 1;
fi;

echo "gathering $vcfDir/*.badsites.txt"
bad="$vcfDir/allsites.txt"
sort $vcfDir/*.badsites.txt | uniq > $bad

echo "vcfToAlignment.pl"
vcfToAlignment.pl bam/*.sorted.bam vcf/*.vcf -o $out -r $ref -b $bad -a 20

echo "convertAlignment.pl"
convertAlignment.pl -i $out -o $out.phy -f phylip -r

