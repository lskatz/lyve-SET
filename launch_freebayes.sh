#!/bin/sh
ref=$1
bam=$2
out_vcf=$3
minAltFrac=$4
minCoverage=$5

readSnpLimit=10;

b=`basename $out_vcf .vcf`;

if [ "$out_vcf" = "" ]; then
  echo usage: $0 ref.fasta query.bam out.vcf min_alt_frac min_coverage;
  echo "  where min_al_frac can be 0.75 and min_coverage can be 10"
  exit 1;
fi;

if [ -e "$out_vcf" ]; then
  echo "Found $out_vcf! I will not rewrite it."
  exit 0;
fi;

# freebayes params
if [ "$minAltFrac" = "" ]; then minAltFrac=0.75; fi;
if [ "$minCoverage" = "" ]; then minCoverage=10; fi;

# for filtering, for later
new_vcf="$out_vcf".tmp
b=`basename $out_vcf .vcf`;
d=`dirname $out_vcf`;
unfiltered_vcf="$d/$b.unfiltered.vcf"
bad=$out_vcf".badsites.txt"

# freebayes
echo "Running Freebayes"
echo "$bam => $out_vcf"

freebayes       \
                `# input and output`\
                --bam $bam \
                --vcf $unfiltered_vcf \
                --fasta-reference $ref \
                \
                `# reporting` \
                --pvar 0.0001 `# Report sites if the probability that there is a polymorphism at the site is greater than N.  default: 0.0001`\
                \
                `# population model`\
                --ploidy 1 `# Sets the default ploidy for the analysis to N.  default: 2`\
                \
                `# allele scope`\
                --no-mnps `# Ignore multi-nuceotide polymorphisms, MNPs.`\
                \
                `# indel realignment`\
                `#--left-align-indels` `# Left-realign and merge gaps embedded in reads. default: false`\
                `#--no-indels` \
                \
                `# input filters`\
                --min-mapping-quality 0 `# Exclude alignments from analysis if they have a mapping quality less than Q.  default: 30`\
                --min-base-quality 20 `# Exclude alleles from analysis if their supporting base quality is less than Q.  default: 20`\
                `#--read-snp-limit $readSnpLimit` `# Exclude reads with more than N base mismatches, ignoring gaps with quality >= mismatch-base-quality-threshold. default: ~unbounded`\
                `# --indel-exclusion-window 5` `# Ignore portions of alignments this many bases from a putative insertion or deletion allele.  default: 0`\
                --min-alternate-fraction $minAltFrac `# Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position.  default: 0.0`\
                --min-coverage $minCoverage `# Require at least this coverage to process a site.  default: 0`

if [ $? -gt 0 ]; then exit 1; fi;

echo "Filtering FreeBayes calls"
filterVcf.pl $unfiltered_vcf --noindels -d $minCoverage -o $new_vcf -b $bad
if [ "$?" -gt 0 ]; then exit 1; fi;
mv -v $new_vcf $out_vcf
if [ "$?" -gt 0 ]; then exit 1; fi;

