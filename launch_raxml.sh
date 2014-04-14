#!/bin/sh
#$ -S /bin/sh
#$ -pe smp 8
#$ -cwd
#$ -V
#$ -o launch_raxml.sh.out -j y
# converts alignment to phylip and then makes a tree out of it

aln=$1
prefix=$2
if [ "$prefix" = "" ]; then
  echo Usage: `basename $0` aln.phy outprefix
  exit 1;
fi

# get the extension
b=`basename $aln`;
suffix="${b##*.}";
if [ "$suffix" != "phy" ]; then
  echo "Converting to phylip format because I did not see a phy extension";
  convertAlignment.pl -f phylip -i $aln -o $aln.phy;
  if [ $? -gt 0 ]; then exit 1; fi;
  aln="$aln.phy";
fi;

raxmlHPC-PTHREADS -f a -s $aln -n $prefix -T 8 -m GTRGAMMA -p $RANDOM -x $RANDOM -N 100
if [ $? -gt 0 ]; then exit 1; fi;

