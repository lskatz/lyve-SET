#!/bin/bash

correctRegion="NC_001416.1:1-48502"
region=$(makeRegions.pl lambda/reference/reference.fasta)
files_to_test="lambda/reference/reference.fasta lambda/bam/sample1.fastq.gz-reference.sorted.bam lambda/vcf/sample2.fastq.gz-reference.vcf.gz"

for file in $files_to_test; do
  region=$(makeRegions.pl $file)
  if [ "$region" != "$correctRegion" ]; then
    echo "ERROR: did not get correct region from $file"
    exit 1
  fi
done
