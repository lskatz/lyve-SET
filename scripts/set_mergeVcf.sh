#!/bin/bash
# Uses bcftools to merge vcf files

script=$(basename $0);
bindir=$(dirname $0);
export PATH=$bindir:$PATH; # enforce the Lyve-SET path

function logmsg(){
  echo "$script: $@"
}

NUMCPUS=1; # default number of cpus
while getopts "o:n:t:r:s" o; do
  case "${o}" in
    o)
      OUT="$OPTARG"
      ;;
    n)
      NUMCPUS="$OPTARG"
      ;;
    t)
      TEMPDIR="$OPTARG"
      ;;
    s)
      ALSOSNPS=1
      ;;
    *)
      logmsg "ERROR: I do not understand option --$o $OPTARG"
      exit 1
      ;;
  esac
done
shift $(($OPTIND-1)) # remove the flag arguments from ARGV

IN="$@";
if [ "$IN" == "" ] || [ "$OUT" == "" ]; then
  echo "$script: merges VCF files that are compressed with bgzip and indexed with tabix"
  echo "Usage: $script -o pooled.vcf.gz *.vcf.gz"
  echo "  -o output.vcf.gz         The output compressed VCF file"
  echo "  -t /tmp/mergeVcfXXXXXX   A temporary directory that will be completely removed upon exit"
  echo "  -s                       Output SNPs only file in addition to the whole thing"
  echo "  -n 1                     Number of cpus"
  exit 1;
fi

SNPSOUT=$(dirname $OUT)/$(basename $OUT .vcf.gz)".snps.vcf.gz"

if [ "$TEMPDIR" == "" ]; then
  TEMPDIR=$(mktemp -d /tmp/mergeVcf.XXXXXX);
  trap "rm -rvf $TEMPDIR" EXIT TERM INT
fi
mkdir -pv "$TEMPDIR"; # just in case
logmsg "temporary directory is $TEMPDIR";

# Run makeRegions.pl based on how many cpus there are 
# and really parallelize this script!
logmsg "Dividing the genome into regions, making it easier to multithread"
makeRegions.pl --numchunks $NUMCPUS $IN > $TEMPDIR/regions.txt
if [ $? -gt 0 ]; then 
  logmsg "ERROR running makeRegions.pl!";
  exit 1;
fi
REGION=$(cat $TEMPDIR/regions.txt);

# Multithread, one thread per region
logmsg "Running bcftools merge";
export IN;
export script;
# TODO if we upgrade past v1.3.1, add in the parameter '--filter-logic x'
echo "$REGION" | xargs -P $NUMCPUS -n 1 bash -c 'echo "$script: merging SNPs in $0"; out='$TEMPDIR'/merged.$$.vcf; bcftools merge --merge all --regions "$0" --force-samples -o $out $IN && bgzip $out && tabix $out.gz && echo "$script: finished with region $0";'
if [ $? -gt 0 ]; then
  logmsg "ERROR with bcftools merge"
  exit 1;
fi

logmsg "Concatenating vcf output"
bcftools concat --allow-overlaps --remove-duplicates $TEMPDIR/merged.*.vcf.gz > $TEMPDIR/concat.vcf
if [ $? -gt 0 ]; then
  logmsg "ERROR with bcftools concat"
  exit 1;
fi

# Generate a SNPs-only merged file, if requested
if [ "$ALSOSNPS" == 1 ]; then
  logmsg "parsing $TEMPDIR/concat.vcf for SNPs"
  bcftools view --include '%TYPE="snp"' $TEMPDIR/concat.vcf > $TEMPDIR/hqPos.vcf
  if [ $? -gt 0 ]; then
    logmsg "ERROR with bcftools view"
    exit 1;
  fi
else
  logmsg "user did not request SNPs-only file also"
fi

# BGzip, tabix, and mv hqSNPs and concat vcfs
for VCF in $TEMPDIR/concat.vcf $TEMPDIR/hqPos.vcf; do
  if [ ! -e "$VCF" ]; then
    continue;
  fi;

  bgzip $VCF
  if [ $? -gt 0 ]; then
    logmsg "ERROR with bgzip $VCF"
    exit 1;
  fi
  tabix $VCF.gz
  if [ $? -gt 0 ]; then
    logmsg "ERROR with tabix $VCF.gz"
    exit 1;
  fi
done

# mv over the SNPs-only files
if [ "$ALSOSNPS" == 1 ]; then
  mv -v $TEMPDIR/hqPos.vcf.gz $SNPSOUT
  mv -v $TEMPDIR/hqPos.vcf.gz.tbi $SNPSOUT.tbi
  logmsg "SNPs-only file can be found in $SNPSOUT";
fi

# Lastly, mv over the large file
mv -v $TEMPDIR/concat.vcf.gz $OUT
mv -v $TEMPDIR/concat.vcf.gz.tbi $OUT.tbi
logmsg "Output file can be found in $OUT"

exit 0;

