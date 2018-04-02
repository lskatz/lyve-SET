#!/bin/bash

set_test.pl --numcpus 2 lambda lambda -- --noqsub
min_snps=$(sort -k3,3n lambda/msa/out.pairwise.tsv | head -n 1 | cut -f 3)
max_snps=$(sort -k3,3n lambda/msa/out.pairwise.tsv | tail -n 1 | cut -f 3)

echo "testing lambda/msa/out.pairwise.tsv"
cat lambda/msa/out.pairwise.tsv

if [ "$min_snps" != "37" ]; then
  echo "I expected the minimum number of SNPs to be 37, but it is $min_snps";
  exit 1
fi
if [ "$max_snps" != "92" ]; then
  echo "I expected the maximum number of SNPs to be 92, but it is $max_snps";
  exit 1
fi

echo "Everything passed!"

