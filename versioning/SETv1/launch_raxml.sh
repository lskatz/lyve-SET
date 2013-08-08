#!/bin/sh

# converts alignment to phylip and then makes a tree out of it

aln=$1
prefix=$2
if [ "$prefix" = "" ]; then
  echo Usage: `basename $0` aln.phy outprefix
  exit 1;
fi

donefile="../progress/raxml.done"

raxmlHPC-PTHREADS -f a -s $aln -n $prefix -T 12 -m GTRGAMMA -p $RANDOM -x $RANDOM -N 100

touch $donefile
