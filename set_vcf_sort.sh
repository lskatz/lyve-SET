#!/bin/bash

# sorts a vcf file
# Author: Lee Katz

VCF=$1;

script=$(basename $0);
if [ "$VCF" = "" ]; then
  echo "$0: sort an uncompressed vcf file";
  echo "Usage: $script in.vcf > sorted.vcf";
  exit 1;
fi

# Get all the headers first
grep '^#' "$VCF"

# For all the lines that aren't headers,
# Sort them by the contig and position
grep -v '^#' "$VCF" | sort -k1,1 -k2,2n 
