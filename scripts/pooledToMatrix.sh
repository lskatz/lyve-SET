#!/bin/bash
# Runs bcftools query to make a SNP matrix

script=$(basename $0);

usage () {
  echo "$script: generates a SNP matrix using BCFtools and a pooled VCF file"
  echo "USAGE: $script -o bcfmatrix.tsv pooled.vcf.gz"
  return 0;
}

logmsg () {
  echo "$script: $@";
  return 0;
}

while getopts "ho:" o; do
  case "${o}" in
    h)
      usage
      exit 1
      ;;
    o)
      OUT="$OPTARG"
      ;;
    *)
      echo "ERROR: I do not understand option $OPTARG"
      ;;
  esac
done
shift $(($OPTIND-1)) # remove the flag arguments from ARGV

IN=$1

if [ "$OUT" == "" ] || [ "$IN" == "" ]; then
  usage
  exit 1;
fi

# Create the basic matrix
if [ ! -e "$IN.tbi" ]; then
  logmsg "$IN.tbi not found -- running tabix."
  tabix $IN
fi

# bcftools query
t='\t'
command="bcftools query -i '%TYPE=\"snp\"' -f '%CHROM$t%POS$t%REF$t[%TGT$t]\\n' --print-header $IN > $OUT.unrefined.tmp"
logmsg $command;
eval $command
if [ $? -gt 0 ]; then
  echo -e "ERROR with bcftools:\n  $command";
  exit 1
fi

# post-process the matrix
logmsg "Post-processing"
perl -lane '
  BEGIN{
    # print the header
    $line=<>;
    chomp($line);
    print $line;
    $numFields=scalar(split(/\t/,$line));
    $lastIndex=$numFields-1;
  }
  for(@F[3..$lastIndex]){
    $_=substr($_,0,1);
  }
  print join("\t",@F);
  ' < $OUT.unrefined.tmp > $OUT.tmp;

mv -v $OUT.tmp $OUT;
rm -v $OUT.unrefined.tmp;
