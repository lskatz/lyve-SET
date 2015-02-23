#!/bin/bash
# Uses bcftools to merge vcf files

script=$(basename $0);

while getopts "o:" o; do
  case "${o}" in
    o)
      OUT="$OPTARG"
      ;;
    *)
      echo "ERROR: I do not understand option $OPTARG"
      ;;
  esac
done
shift $(($OPTIND-1)) # remove the flag arguments from ARGV

IN="$@";
if [ "$IN" == "" ] || [ "$OUT" == "" ]; then
  echo "$script: merges VCF files that are compressed with bgzip and indexed with tabix"
  echo "Usage: $script -o pooled.vcf.gz *.vcf.gz"
  echo "  -o output.vcf.gz  The output compressed VCF file"
  exit 1;
fi

command="bcftools merge $IN -O z > $OUT.tmp && mv -v $OUT.tmp $OUT"
eval $command
if [ $? -gt 0 ]; then 
  echo -e "ERROR with bcftools:\n  $command";
  rm -vf $OUT.tmp
  exit 1;
fi;

tabix $OUT
if [ $? -gt 0 ]; then
  echo "ERROR with tabix"
  exit 1;
fi;

