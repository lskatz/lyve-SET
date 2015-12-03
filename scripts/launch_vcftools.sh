#!/bin/bash


set -e
function usage() { 
	echo "
	$0: takes an mpileup+ref or vcf file and uses VCFtools and other VCFlib binaries to produce a vcf file with hqSNPs
	Usage: `basename $0` [-m file.mpileup && -r ref.fasta] || [-v file.vcf[.gz]]
		--tempdir     ./tmp/  A temporary directory to store files
		--coverage         4  Min coverage
		--allowedFlanking  0  Neighboring bases to filter between SNPs
		--altfreq       0.75  Min consensus agreement [%] for a SNP
		--minQ            30  Minimum average Phred value for a SNP
		--region    file.bed  File of positions to include (default without bed file: read all positions)
		--exclude   file.bed  File of positions to mask    (conflicts with --region; default: no site exclusions)

		--xopts '' [optional] Extra options to pass to VCFtools
		
	Paths can be absolute or relative.
	"
	}


# defaults
ALTFREQ='--non-ref-af 0.75'
COVERAGE='--minDP 4'
MINQ='--minQ 30'
TMP=tmp
VCFFMT='--vcf'


nopts=$#
for ((i=1; i <= nopts; i++)); do
	case "$1" in
		-a | --allowedFlanking) # adjacent SNPs filtered out
			FLANK="--thin $2"
			echo "    Discarding SNP sites $2 bases apart"
			shift 2
			;;
		-c | --coverage)  # minimum coverage depth
			COVERAGE="--minDP $2"
			echo "    Min coverage: $2"
			shift 2
			;;
		-e | --exclude)  # BED coordinates to mask
			EXCLUDE="--exclude-bed $2"
			echo "    Excluding coordinates from: $2"
			shift 2
			;;
		-f | --altfreq)  # minimum read consensus [%] for a SNP
			ALTFREQ="--non-ref-af $2"
			echo "    Min consensus: $2"
			shift 2
			;;
		-i | --region)  # BED coordinates to include
			REGION="--bed $2"
			echo "    Including coordinates from: $2"
			shift 2
			;;
		-m | --mpileup)  # input mpileup file
			MPILEUP="$2"
			echo "    Input mpileup file: $2"
			shift 2
			;;
		-n | --numcpus)  # VCFtools not threaded
			NUMCPUS="$2"
			echo "    CPUs requested: $2"
			shift 2
			;;
		-o | --outpref)
			OUTPREF="--out $2"
			echo "    Output prefix: $2"
			shift 2
			;;
		-r | --reference)
			REF="$2"
			echo "    Genome reference: $2"
			shift 2
			;;
		-t | --tempdir)
			TMP="$2"
			echo "    Temporary directory (path): $2"
			shift 2
			;;
		-v | --vcf)
			VCF="$2"
			echo "    VCF file: $2"
			shift 2
			;;
		-q | --minQ)
			MINQ="--minQ $2"
			echo "    Minimum average Phred score: $2"
			shift 2
			;;
		-x | --xopts)
			XOPTS="$2"
			echo "    extra options passed to VCFtools: $2"
			shift 2
			;;
		-h | --help)
			usage
			exit 1
			;;
	esac
done


# Input file and dependency checks
[[ -z "$VCF" ]] && [[ -z "$MPILEUP" ]] && { usage; exit 1; }
[[ -n "$MPILEUP" ]] && [[ -z "$REF" ]] && { echo 'ERROR: mpileup requires reference genome (-r) input as well'; usage; exit 1; }
[[ -n "$REF" ]] && [[ -z "$MPILEUP" ]] && { echo 'ERROR: reference genome requires mpileup (-m) input as well'; usage; exit 1; }
command -v vcffilter >/dev/null 2>&1 || { echo 'ERROR: vcffilter binary not found' >&2; exit 1; }
command -v vcffixup >/dev/null 2>&1 || { echo 'ERROR: vcffixup binary not found' >&2; exit 1; }
command -v vcf-sort >/dev/null 2>&1 || { echo 'ERROR: vcf-sort script not found' >&2; exit 1; }
command -v vcftools >/dev/null 2>&1 || { echo 'ERROR: vcftools binary not found' >&2; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo 'ERROR: bgzip binary not found' >&2; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo 'ERROR: tabix binary not found' >&2; exit 1; }

# Special exception when both '--region' and '--exclude' are requested
if [[ -n "$REGION" ]] && [[ -n "$EXCLUDE" ]]; then
	MINUEND=$(echo "$REGION" | sed 's/^--region //g')
	SUBTRAHEND=$(echo "$EXCLUDE" | sed 's/^--exclude-bed //g')
	U=$(echo "$EXCLUDE" | sed -r 's/(\.bed$|^--exclude-bed )//g')
	#substract the two conflicting opts --exclude --region' for vcftools
	command -v bedtools >/dev/null 2>&1 || { echo 'ERROR: bedtools binary not found' >&2; exit 1; }
	bedtools subtract -a "$MINUEND" -b "$SUBTRAHEND" > "$TMP"/"$U".difference.bed
	DIFF="--bed $TMP/$U.difference.bed"
	REGION=''
	EXCLUDE=''
fi

# Manage tmp dir
if [ ! -d "$TMP" ]; then
	mkdir -p "$TMP"
	echo '    Created temporary directory...'
else
	echo "WARNING: temp dir $TMP already exists and conflicting files will be overwritten"
fi

# convert mpileup to vcf when mpileup file given and vcf not given
if [[ -n "$MPILEUP" ]]; then
	command -v bcftools >/dev/null 2>&1 || { echo 'ERROR: bcftools binary not found' >&2; exit 1; }
	command -v samtools >/dev/null 2>&1 || { echo 'ERROR: samtools binary not found' >&2; exit 1; }
	b=$(basename "$MPILEUP" .sorted.bam)
	samtools faidx "$REF"
	#samtools view -H "$MPILEUP" | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P "$NUMCPUS" sh -c "samtools mpileup -d 1000 -suvf $REF -r '{}' $MPILEUP 1> $TMP/$b.vcf 2> /dev/null" 
	samtools mpileup -d 1000 -suvf "$REF" "$MPILEUP" 2> /dev/null | bcftools call -c > "$TMP"/"$b".vcf
	VCF="$TMP"/"$b".vcf
	echo 'mpileup conversion into VCF completed'
fi

if [[ -z "$OUTPREF" ]]; then
	if [[ -n "$VCF" ]]; then
		B=$(basename "$VCF" | sed -r 's/\.(vcf|vcf\.gz)$//1')
		OUTPREF="--out ${TMP}/${B}"
	else
		B=$(basename "$MPILEUP" | sed -r 's/\.(bam|sorted\.bam|mpileup)$//1')
		OUTPREF="--out ${TMP}/${B}"
	fi
fi

[[ "$VCF" == *.gz ]] && VCFFMT='--gzvcf' #enables compressed or uncompressed

# run VCFtools
file=${VCF##*/}
b=${file%%.*}
echo "$b" > "$TMP"/sampleID
COMMAND="$(echo $VCFFMT $VCF $OUTPREF $COVERAGE $FLANK $MINQ $REGION $EXCLUDE $DIFF \
    --remove-indels --remove-filtered-all $XOPTS --recode --recode-INFO-all \
    | sed 's/ \{1,\}/ /g')" #cleanup spaces
vcftools $COMMAND 2> /dev/null
echo 'VCFtools completed'
[[ -n "$MPILEUP" ]] && rm -v "$TMP"/"$B".vcf
#discard PL field due to inconsisent lines which breaks downstream mergeVcf.sh (that invokes bcftools merge); some lines='GT:PL	1/1:255,63,0', others='GT:PL	1/1:0,105:94:99:255,255,0'
vcffilter -g 'GT = 1/1' "$TMP"/"$B".recode.vcf | vcffixup - | vcffilter -f 'AC > 0' | vcf-sort -c | bcftools annotate -x PL,FORMAT - | bcftools reheader -s "$TMP"/sampleID - > "$TMP"/"$B".vcf 
echo 'vcffilter completed'

if grep -Pq '^#CHROM[A-Z\t]+sm$' "$TMP"/"$B".vcf; then
	#clean up sampleID in VCF file
	b=$(basename "$B" .fastq-reference)
	sed -r -i "/^#CHROM[A-Z[[:space:]]+sm$/s/sm/$b/1" "$TMP"/"$B".vcf
	echo "changed sample ID inside VCF to reflect input filename: $b"
fi

bgzip -cf "$TMP"/"$B".vcf > "$TMP"/"$B".vcf.gz
rm -v "$TMP"/"$B".recode.vcf "$TMP"/"$B".vcf "$TMP"/sampleID
[[ -z "$DIFF" ]] && { rm -v "$TMP"/"$U".difference.bed }
tabix -f "$TMP"/"$B".vcf.gz #&> /dev/null
echo 'bgzip compression and tabix indexing completed'
#to-do: filter for at least 1 read on both -f "ADF > 0 & ADR > 0" doesn't work

