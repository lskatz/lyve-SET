#!/bin/bash

script=$(basename $0)

## Flag arguments
numcpus=1
outdir="out"
while getopts ":n:o:" o; do
  case "${o}" in
    n)
      numcpus=$OPTARG
      ;;
    o)
      outdir=$OPTARG
      ;;
    *)
      echo "ERROR: I do not understand option $OPTARG"
      ;;
  esac
done
shift $(($OPTIND-1)) # remove the flag arguments from ARGV

if [ "$1" == "" ]; then
  echo "Shuffles a set of reads into an output directory"
  echo "Usage: $0 -o outdir *.fastq.gz"
  echo "  -n numcpus";
  exit 1;
fi;

# Pair up everything
ARGS=""
for f ; do
  if [[ "$f" =~ R2 ]]; then
    continue;
  fi;
  b=$(sed 's/_R1_.*$//' <<< $f);
  r=$(ls ${b}*_R2_*.fastq.gz)
  ARGS="$f $r $outdir/$b.fastq.gz $ARGS";
done;
ARGS=$(sed 's/^ \| $//g' <<< $ARGS); # whitespace trimming just in case

# Multithread the shuffling
mkdir -p $outdir
echo "$ARGS" | \
  xargs -P "$numcpus" -n 3 sh -c '
    echo "Shuffling into $2";
    if [ -e "$2" ]; then 
      echo Found $2
      exit 0; 
    fi; 
    run_assembly_shuffleReads.pl $0 $1 | gzip -c > $2.tmp && mv $2.tmp $2
    echo "Finished shuffling into $2"
  '

exit 0;
