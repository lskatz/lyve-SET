#!/bin/bash
# Uses bcftools to merge vcf files

script=$(basename $0);

NUMCPUS=1; # default number of cpus
while getopts "o:n:" o; do
  case "${o}" in
    o)
      OUT="$OPTARG"
      ;;
    n)
      NUMCPUS="$OPTARG"
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

# Figure out regions names using the VCF index files.
# Do not multithread xargs here because of race conditions and stdout.
REGION=$(echo "$IN" | xargs -P 1 -n 1 tabix -l | sort | uniq);

# Multithread, one thread per region
TMPDIR=$(mktemp -d /tmp/mergeVcf.XXXXXX);
echo "$script: temporary directory is $TMPDIR";
echo "$script: Running bcftools merge";
export IN;
echo "$REGION" | xargs -P $NUMCPUS -n 1 -I {} bash -c 'bcftools merge --regions "{}" $IN --force-samples -o '$TMPDIR'/merged.$$.vcf;'
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bcftools merge"
  rm -rvf $TMPDIR;
  exit 1;
fi

echo "$script: Concatenating vcf output"
bcftools concat $TMPDIR/merged.*.vcf > $TMPDIR/concat.vcf
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bcftools concat"
  rm -rvf $TMPDIR;
  exit 1;
fi

echo "$script: parsing $TMPDIR/concat.vcf for high-quality positions"
perl -lane '
    if(/^#/){
      print;
    }elsif($F[4] ne "N"){ 
      print; 
    }else { 
      $numF=@F; 
      $has_hq=0; 
      for my $info(@F[9..$numF-1]){ 
        $gt=substr((split(/:/,$info))[0],0,1); 
        $has_hq=1 if($gt=~/[ATCG0\.,]/i);
      } 
      print if($has_hq); 
    }
' < $TMPDIR/concat.vcf > $TMPDIR/hqPos.vcf
if [ $? -gt 0 ]; then
  echo "$script: ERROR with perl parsing"
  rm -rvf $TMPDIR;
  exit 1;
fi

bgzip $TMPDIR/hqPos.vcf
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bgzip"
  rm -rvf $TMPDIR;
  exit 1;
fi
tabix $TMPDIR/hqPos.vcf.gz
if [ $? -gt 0 ]; then
  echo "$script: ERROR with tabix"
  rm -rvf $TMPDIR;
  exit 1;
fi

mv -v $TMPDIR/hqPos.vcf.gz $OUT
mv -v $TMPDIR/hqPos.vcf.gz.tbi $OUT.tbi

rm -rvf $TMPDIR;

exit 0;


# zcat out.pooled.vcf.gz|perl -lane 'if(/^#/){print;} elsif($F[4] ne "N"){ print; } else { $numF=@F; $has_hq=0; for my $info(@F[8..$numF-1]){ $gt=substr((split(/:/,$info))[0],0,1); $has_hq=1 if($gt=~/[ATCG0\.,]/i);} print if($has_hq); }' | grep -cv '#'

# This command does the following:
# 1. bcftools merge (merging VCF samples into one)
# 2. Checks to see if the site is high-quality or not
#   a. Is it a header or is there no ambiguous alt call? -- print.
#   b. Is there at least one sample with an unambiguous base call? --print.
#   c. Are all bases ambiguous? -- skip
#command="bcftools merge $IN --force-samples | \
#perl -lane 'if(/^#/){print;} elsif(\$F[4] ne \"N\"){ print; } else { \$numF=@F; \$has_hq=0; for my \$info(@F[9..\$numF-1]){ \$gt=substr((split(/:/,\$info))[0],0,1); \$has_hq=1 if(\$gt=~/[ATCG0\.,]/i);} print if(\$has_hq); }' |\
#bgzip -c > $OUT.tmp
#"
#eval $command
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

