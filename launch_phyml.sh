#!/bin/sh

# makes a tree out of an aln

aln=$1
if [ "$aln" = "" ]; then
  echo Usage: `basename $0` aln.phy
  exit 1;
fi

# find which phyml to use
phyml=`(which phyml || which PhyML || which PhyML-3.1_linux64 || which phyml_linux_64) 2>/dev/null`;
if [ $? -gt 0 ]; then echo "Could not find phyml"; exit 1; fi;
echo "Found phyml at $phyml"

$phyml -i $aln -b -4 -m GTR -s BEST --quiet
if [ $? -gt 0 ]; then exit 1; fi;

