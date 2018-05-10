#!/bin/bash


min_snps=$(sort -k3,3n lambda/msa/out.pairwise.tsv | head -n 1 | cut -f 3)
max_snps=$(sort -k3,3n lambda/msa/out.pairwise.tsv | tail -n 1 | cut -f 3)

if [ "$min_snps" != "37" ]; then
  echo "I expected the minimum number of SNPs to be 37, but it is $min_snps";
  cat lambda/msa/out.pairwise.tsv
  exit 1
fi
if [ "$max_snps" != "92" ]; then
  echo "I expected the maximum number of SNPs to be 92, but it is $max_snps";
  cat lambda/msa/out.pairwise.tsv
  exit 1
fi


