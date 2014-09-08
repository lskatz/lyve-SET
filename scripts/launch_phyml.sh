#!/bin/sh
#$ -S /bin/sh
#$ -pe smp 1
#$ -cwd
#$ -V
#$ -o launch_phyml.sh.out -j y

# makes a tree out of an aln

script=`basename $0`;
aln=$1
if [ "$aln" = "" ]; then
  echo Usage: `basename $0` aln.phy
  exit 1;
fi

# find which phyml to use
phyml=`(which phyml || which PhyML || which PhyML-3.1_linux64 || which phyml_linux_64) 2>/dev/null`;
if [ $? -gt 0 ]; then echo "Could not find phyml"; exit 1; fi;
echo "$script: Found phyml at $phyml"

# get the extension
b=`basename $aln`;
suffix="${b##*.}";
if [ "$suffix" != "phy" ]; then
  echo "Converting to phylip format because I did not see a phy extension";
  convertAlignment.pl -f phylip -i $aln -o $aln.phy;
  if [ $? -gt 0 ]; then exit 1; fi;
  aln="$aln.phy";
fi;


yes | \
  $phyml -i $aln -b -4 -m GTR -s BEST --quiet;
if [ $? -gt 0 ]; then echo "$script: ERROR in phyml" exit 1; fi;
echo "$script: Finished without error!";

