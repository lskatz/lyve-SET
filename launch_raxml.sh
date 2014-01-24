#!/bin/sh
#$ -V 
#$ -cwd
#$ -o raxml.log
#$ -j y
#$ -pe smp 8

# converts alignment to phylip and then makes a tree out of it

aln=$1
prefix=$2
if [ "$prefix" = "" ]; then
  echo Usage: `basename $0` aln.phy outprefix
  exit 1;
fi

raxmlHPC-PTHREADS -f a -s $aln -n $prefix -T 8 -m GTRGAMMA -p $RANDOM -x $RANDOM -N 100
if [ $? -gt 0 ]; then exit 1; fi;

