#!/bin/bash


set -e
function usage() { 
	echo "
	Usage: `basename $0` [-m file.mpileup && -r ref.fasta] || [-v file.vcf[.gz]]
		--tempdir     ./tmp/  A temporary directory to store files
		--coverage        10  Min coverage
		--allowedFlanking  0  Neighboring bases to filter between SNPs
		--altfreq       0.75  Min consensus agreement [%] for a SNP
		--minQ            30  Minimum average Phred value for a SNP
		--region    file.bed  File of positions to include (default without bed file: read all positions)
		--exclude   file.bed  File of positions to mask    (conflicts with --region; default: no site exclusions)

		--xopts '' [optional] Extra options to pass to vcftools
		
	Paths can be absolute or relative.
	"
	}


# defaults
ALTFREQ='--non-ref-af 0.75'
COVERAGE='--minDP 10'
#MINQ='--minQ 30'
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
		-n | --numcpus)  # vcftools not threaded
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
		# -q | --minQ)
		# 	MINQ="--minQ $2"
		# 	echo "    Minimum average Phred score: $2"
		# 	shift 2
		# 	;;
		-x | --xopts)
			XOPTS="$2"
			echo "    extra options passed to vcftools: $2"
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
[[ -n "$REGION" ]] && [[ -n "$EXCLUDE" ]] && { echo 'ERROR: conflicting opts --exclude --region'; usage; exit 1; }
command -v vcftools >/dev/null 2>&1 || { echo 'ERROR: vcftools binary not found' >&2; exit 1; }
command -v vcffilter >/dev/null 2>&1 || { echo 'ERROR: vcffilter binary not found' >&2; exit 1; }

# Manage tmp dir
if [ ! -d "$TMP" ]; then
	mkdir -p "$TMP"
	echo '    Created temporary directory...'
else
	echo "WARNING: temp dir $TMP already exists and conflicting files will be overwritten"
fi

# convert mpileup to vcf
if [[ -n "$MPILEUP" ]]; then
	b=$(basename "$MPILEUP" .sorted.bam)
	samtools faidx "$REF"
	samtools view -H "$MPILEUP" | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P "$NUMCPUS" sh -c "samtools mpileup -d 1000 -suvf $REF -r '{}' $MPILEUP 1> $TMP/$b.vcf 2> /dev/null"
	#samtools mpileup -d 1000 -suvf "$REF" "$MPILEUP" 2> /dev/null > "$TMP"/"$b".vcf
	#tabix "$TMP"/"$b".vcf.gz &> /dev/null
	VCF="$TMP"/"$b".vcf
fi

if [[ -z "$OUTPREF" ]]; then
	B=$(basename "$VCF" | sed -r 's/\.(vcf|vcf\.gz)$//1')
	OUTPREF="--out ${TMP}/${B}"
fi

[[ "$VCF" == *.gz ]] && VCFFMT='--gzvcf' #enables compressed or uncompressed

# run vcftools
B=$(basename "$VCF" | sed -r 's/\.(vcf|vcf\.gz)$//1')
COMMAND="$(echo $VCFFMT $VCF $OUTPREF $COVERAGE $FLANK $MINQ $REGION $EXCLUDE \
    --remove-indels --remove-filtered-all $XOPTS --recode --recode-INFO-all \
    | sed 's/ \{1,\}/ /g')" #cleanup spaces
echo 'mpileup complete'
vcftools "$COMMAND" &> /dev/null
echo 'vcftools complete'
#rm -v "$TMP"/"$B".log
vcffilter -g "GT = 1/1" "$TMP"/"$B".recode.vcf | vcffixup - | vcffilter -f "AC > 0" | bgzip -f > "$TMP"/"$B".vcf.gz
echo 'vcffilter complete'
tabix -f "$TMP"/"$B".vcf.gz #&> /dev/null
echo 'tabix complete'
#rm -v "$TMP"/"$B".recode.vcf
#to-do: filter for at least 1 read on both -f "ADF > 0 & ADR > 0" doesn't seem to work

