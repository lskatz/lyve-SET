#!/bin/bash
# converts alignment to phylip and then makes a tree out of it
# Author: Lee Katz <lkatz@cdc.gov>
# TODO: make an outgroup parameter, probably using eigenvectors (least connected is the best guess for an outgroup)

#$ -S /bin/bash
#$ -pe smp 2
#$ -cwd
#$ -V
#$ -o launch_raxml.sh.out -j y

aln=$1
prefix=$2
if [ "$prefix" = "" ]; then
  echo Usage: `basename $0` aln.phy outprefix
  exit 1;
fi

# If the user specifies the outgroup, then set the parameter for Raxml
outgroup=$3
outgroupParam=""
if [ "$outgroup" = "" ]; then
  echo Warning: outgroup was not given.
else
  outgroupParam="-o $outgroup"
fi

# set the number of CPUs according to the number of slots, or 8 if undefined
numcpus=${NSLOTS:-2}

# get the extension
b=`basename $aln`;
suffix="${b##*.}";
if [ "$suffix" != "phy" ]; then
  echo "Converting to phylip format because I did not see a phy extension";
  convertAlignment.pl -f phylip -i $aln -o $aln.phy;
  if [ $? -gt 0 ]; then exit 1; fi;
  aln="$aln.phy";
fi;

# Which raxml to use?
which raxmlHPC 1>/dev/null 2>&1 && EXEC=raxmlHPC
which raxmlHPC-PTHREADS 1>/dev/null 2>&1 && EXEC=raxmlHPC-PTHREADS

# what version is raxml?
VERSION=$($EXEC -v | grep -o 'version [0-9]' | grep -o [0-9]);
echo "I detect that you are using version $VERSION of raxml."

if [ $VERSION = 8 ]; then
  MODEL=ASC_GTRGAMMA
else
  MODEL=GTRGAMMA
fi
$EXEC -f a -s $aln -n $prefix -T $numcpus -p $RANDOM -x $RANDOM -N 100 -m $MODEL $outgroupParam
if [ $? -gt 0 ]; then exit 1; fi;

