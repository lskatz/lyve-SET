#!/bin/bash
# Uses bcftools to merge vcf files

script=$(basename $0);

NUMCPUS=1; # default number of cpus
while getopts "o:n:t:" o; do
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
  echo "  -t tempdir        A temporary directory (Default: a random dir using mktemp)"
  exit 1;
fi

# Figure out regions names using the VCF index files.
# Do not multithread xargs here because of race conditions and stdout.
REGION=$(echo "$IN" | xargs -P 1 -n 1 tabix -l | sort | uniq);

# Multithread, one thread per region
if [ "$TEMPDIR" == "" ]; then
  TEMPDIR=$(mktemp -d /tmp/mergeVcf.XXXXXX);
fi
echo "$script: temporary directory is $TEMPDIR";
echo "$script: Running bcftools merge";
export IN;
echo "$REGION" | xargs -P $NUMCPUS -n 1 -I {} bash -c 'bcftools merge --regions "{}" $IN --force-samples -o '$TEMPDIR'/merged.$$.vcf;'
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bcftools merge"
  rm -rvf $TEMPDIR;
  exit 1;
fi

echo "$script: Concatenating vcf output"
bcftools concat $TEMPDIR/merged.*.vcf > $TEMPDIR/concat.vcf
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bcftools concat"
  rm -rvf $TEMPDIR;
  exit 1;
fi

#echo "$script: parsing $TEMPDIR/concat.vcf for high-quality positions"
#perl -lane '
#    if(/^#/){
#      print;
#    }elsif($F[4] ne "N"){ 
#      print; 
#    }else { 
#      $numF=@F; 
#      $has_hq=0; 
#      for my $info(@F[9..$numF-1]){ 
#        $gt=substr((split(/:/,$info))[0],0,1); 
#        $has_hq=1 if($gt=~/[ATCG0\.,]/i);
#      } 
#      print if($has_hq); 
#    }
#' < $TEMPDIR/concat.vcf > $TEMPDIR/hqPos.vcf
#if [ $? -gt 0 ]; then
#  echo "$script: ERROR with perl parsing"
#  rm -rvf $TEMPDIR;
#  exit 1;
#fi

echo "$script: parsing $TEMPDIR/concat.vcf for SNPs"
bcftools annotate --include '%TYPE="snp"' < $TEMPDIR/concat.vcf > $TEMPDIR/hqPos.vcf
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bcftools annotate"
  rm -rvf $TEMPDIR;
  exit 1;
fi

bgzip $TEMPDIR/hqPos.vcf
if [ $? -gt 0 ]; then
  echo "$script: ERROR with bgzip"
  rm -rvf $TEMPDIR;
  exit 1;
fi
tabix $TEMPDIR/hqPos.vcf.gz
if [ $? -gt 0 ]; then
  echo "$script: ERROR with tabix"
  rm -rvf $TEMPDIR;
  exit 1;
fi

mv -v $TEMPDIR/hqPos.vcf.gz $OUT
mv -v $TEMPDIR/hqPos.vcf.gz.tbi $OUT.tbi

rm -rvf $TEMPDIR;

exit 0;

