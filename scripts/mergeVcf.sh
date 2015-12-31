#!/bin/bash
# Uses bcftools to merge vcf files

script=$(basename $0);
bindir=$(dirname $0);
export PATH=$bindir:$PATH; # enforce the Lyve-SET path

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
      echo "ERROR: I do not understand option $OPTARG"
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
  echo "  -n 1                     Number of CPUs"
  echo "  -s                       Output SNPs only file in addition to the whole thing"
  exit 1;
fi

SNPSOUT=$(dirname $OUT)/$(basename $OUT .vcf.gz)".snps.vcf.gz"

if [ "$TEMPDIR" == "" ]; then
  TEMPDIR=$(mktemp --directory /tmp/mergeVcf.XXXXXX);
fi
mkdir -pv "$TEMPDIR"; # just in case
export TEMPDIR;

REGIONSFILE="$TEMPDIR/regions.txt";
touch $REGIONSFILE; #make sure this file is accessible
if [ $? -gt 0 ]; then echo "ERROR: could not access $REGIONSFILE"; exit 1; fi;

# Run makeRegions.pl based on how many cpus there are 
# and really parallelize this script!
echo "$script: Dividing the genome into regions, making it easier to multithread. Regions will be found in $REGIONSFILE"
MAKEREGIONS="makeRegions.pl --numchunks $NUMCPUS --numcpus $NUMCPUS $IN > $REGIONSFILE"
echo $MAKEREGIONS
eval $MAKEREGIONS
if [ $? -gt 0 ]; then 
  echo "$script: ERROR running makeRegions.pl!";
  rm -rvf $TEMPDIR;
  exit 1;
fi
REGION=$(cat $REGIONSFILE);

# IF makeRegions.pl didn't work or wasn't invoked:
# Figure out regions names using the VCF index files.
# Do not multithread xargs here because of race conditions and stdout.
if [ "$REGION" == "" ]; then 
  REGION=$(echo "$IN" | xargs -P 1 -n 1 tabix -l | sort | uniq);
fi


# Multithread, one thread per region
echo "$script: temporary directory is $TEMPDIR";
echo "$script: Running bcftools merge";
export IN;
export script;
echo "$REGION" | xargs -P $NUMCPUS -n 1 -I {} bash -c 'echo "$script: merging SNPs in {}"; out='$TEMPDIR'/merged.$$.vcf; bcftools merge --merge all --regions "{}" --force-samples -o $out $IN && bgzip $out && tabix $out.gz && echo "$script: finished with region {}";'
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bcftools merge"
  rm -rvf $TEMPDIR;
  exit 1;
fi

echo "$script: Concatenating vcf output"
bcftools concat $TEMPDIR/merged.*.vcf.gz | vcf-sort > $TEMPDIR/concat.vcf
#bcftools concat --allow-overlaps --remove-duplicates $TEMPDIR/merged.*.vcf.gz > $TEMPDIR/concat.vcf
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bcftools concat"
  rm -rvf $TEMPDIR;
  exit 1;
fi

# Generate a SNPs-only merged file, if requested
if [ "$ALSOSNPS" == 1 ]; then
  echo "$script: parsing $TEMPDIR/concat.vcf for SNPs"
  bcftools annotate --include '%TYPE="snp"' $TEMPDIR/concat.vcf > $TEMPDIR/hqPos.vcf
  if [ $? -gt 0 ]; then
    echo "$script: ERROR with bcftools annotate"
    rm -rvf $TEMPDIR;
    exit 1;
  fi
else
  echo "$script: user did not request SNPs-only file also"
fi

# BGzip, tabix, and mv hqSNPs and concat vcfs
for VCF in $TEMPDIR/concat.vcf $TEMPDIR/hqPos.vcf; do
  if [ ! -e "$VCF" ]; then
    continue;
  fi;

  echo "$script: Sorting $VCF and removing indels";
  bcftools annotate --include '%TYPE!="indel"' < $VCF > $VCF.tmp && mv $VCF.tmp $VCF;

  echo "$script: compressing $VCF";
  bgzip $VCF
  if [ $? -gt 0 ]; then
    echo "$script: ERROR with bgzip $VCF"
    rm -rvf $TEMPDIR;
    exit 1;
  fi
  echo "$script: indexing $VCF";
  tabix $VCF.gz
  if [ $? -gt 0 ]; then
    echo "$script: ERROR with tabix $VCF.gz"
    rm -rvf $TEMPDIR;
    exit 1;
  fi
done

# mv over the SNPs-only files
if [ "$ALSOSNPS" == 1 ]; then
  mv -v $TEMPDIR/hqPos.vcf.gz $SNPSOUT
  mv -v $TEMPDIR/hqPos.vcf.gz.tbi $SNPSOUT.tbi
  echo "$script: SNPs-only file can be found in $SNPSOUT";
fi

# Lastly, mv over the large file
mv -v $TEMPDIR/concat.vcf.gz $OUT
mv -v $TEMPDIR/concat.vcf.gz.tbi $OUT.tbi
echo "$script: Output file can be found in $OUT"

rm -rvf $TEMPDIR;

exit 0;
