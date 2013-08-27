#!/bin/sh

bamDir=$1
vcfDir=$2
ref=$3
out=$4

if [ "$out" = "" ]; then
  echo usage: $0 bam/ vcf/ ref.fasta out.aln.fas
  echo "  bam directory's *.sorted.bam will be read"
  echo "  vcf dir's *.vcf will be read"
  echo "  vcfDir/*.badsites.txt will be used as 'bad sites', as produced by filterVcf.pl"
  exit 1;
fi;

echo "gathering $vcfDir/*.badsites.txt"
bad="$vcfDir/allsites.txt"
sort $vcfDir/*.badsites.txt | uniq > $bad
if [ $? -gt 0 ]; then exit 1; fi;

echo "vcfToAlignment.pl"
vcfToAlignment.pl $bamDir/*.sorted.bam $vcfDir/*.vcf -o $out -r $ref -b $bad -a 0 -n 8
if [ $? -gt 0 ]; then exit 1; fi;

echo "convertAlignment.pl"
convertAlignment.pl -i $out -o $out.phy -f phylip -r
if [ $? -gt 0 ]; then exit 1; fi;

