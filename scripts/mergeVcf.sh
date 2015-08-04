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

# Figure out regions names using the VCF index files
REGION=$(echo "$IN" | xargs -n 1 tabix -l | sort | uniq);
#echo "$REGION";
#exit 1;

# TODO: multithread, one thread per region


# zcat out.pooled.vcf.gz|perl -lane 'if(/^#/){print;} elsif($F[4] ne "N"){ print; } else { $numF=@F; $has_hq=0; for my $info(@F[8..$numF-1]){ $gt=substr((split(/:/,$info))[0],0,1); $has_hq=1 if($gt=~/[ATCG0\.,]/i);} print if($has_hq); }' | grep -cv '#'

# This command does the following:
# 1. bcftools merge (merging VCF samples into one)
# 2. Checks to see if the site is high-quality or not
#   a. Is it a header or is there no ambiguous alt call? -- print.
#   b. Is there at least one sample with an unambiguous base call? --print.
#   c. Are all bases ambiguous? -- skip
command="bcftools merge $IN --force-samples | \
perl -lane 'if(/^#/){print;} elsif(\$F[4] ne \"N\"){ print; } else { \$numF=@F; \$has_hq=0; for my \$info(@F[9..\$numF-1]){ \$gt=substr((split(/:/,\$info))[0],0,1); \$has_hq=1 if(\$gt=~/[ATCG0\.,]/i);} print if(\$has_hq); }' |\
bgzip -c > $OUT.tmp
"
eval $command
if [ $? -gt 0 ]; then
  echo -e "ERROR with bcftools:\n  $command";
  rm -vf $OUT.tmp
  exit 1;
fi;
mv -v $OUT.tmp $OUT
if [ $? -gt 0 ]; then
  echo "ERROR moving $OUT.tmp to $OUT";
  exit 1;
fi;

tabix $OUT
if [ $? -gt 0 ]; then
  echo "ERROR with tabix"
  exit 1;
fi;

