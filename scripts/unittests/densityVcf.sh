#!/bin/bash

# This test assumes that set_test.pl has already been run
# and that it generated project directory lambda.

####### density filtering ##########
# All positions
positions=$(zgrep -v '^#' lambda/msa/out.pooled.snps.vcf.gz | cut -f 2);
# Density filter of 2 that removes the first site: no two sites should
# be within 158 positions of each other.
filteredPositionsWindow2=$(densityFilterVcf.pl --window 158 --density-filter 2 lambda/msa/out.pooled.snps.vcf.gz | grep -v '^#' | cut -f 2)

smallestDifference2=$(
  echo "$filteredPositionsWindow2" | perl -e '
    @site=<>; 
    chomp(@site); 
    for($i=0;$i<@site-1;$i++){
      print $site[$i+1] - $site[$i]."\n";
    }
  ' |\
  sort -n | head -n 1
)
if [ "$smallestDifference2" -lt 158 ]; then
  echo "Did not remove sites close together than 158 when density filtering was 2"
  echo "Smallest distance between sites was $smallestDifference2"
  exit 1
fi

# Density filter of 3 that removes the first site: no two sites, separated
# by one site, should be within 184 positions of each other.
filteredPositionsWindow3=$(densityFilterVcf.pl --window 184 --density-filter 3 lambda/msa/out.pooled.snps.vcf.gz | grep -v '^#' | cut -f 2)
smallestDifference3=$(
  echo "$filteredPositionsWindow3" | perl -e '
    @site=<>; 
    chomp(@site); 
    for($i=0;$i<@site-2;$i++){
      print $site[$i+2] - $site[$i]."\n";
    }
  ' |\
  sort -n | head -n 1
)
if [ "$smallestDifference3" -lt 184 ]; then
  echo "Did not remove sites close together than 158 when density filtering was 2"
  echo "Smallest distance between sites was $smallestDifference3"
  exit 1
fi

