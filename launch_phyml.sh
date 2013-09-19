#!/bin/sh

# makes a tree out of an aln

aln=$1
if [ "$aln" = "" ]; then
  echo Usage: `basename $0` aln.phy
  exit 1;
fi

PhyML -i $aln -b -4 -m GTR -s BEST --quiet
if [ $? -gt 0 ]; then exit 1; fi;

