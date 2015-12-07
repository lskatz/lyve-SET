#!/bin/bash


set -e
function usage() { 
	echo "
	Usage: `basename $0` -f file.fastq[.gz] -b file.bam -t tmp/ -r reference.fasta [-x '']
		
		-t [default: ./tmp] Path to temporary directory
		-n [default: 1] Number of threads to use
		-pe <0|1|2> [default: 0] For 'auto', single-end, or paired-end respectively
		
		-x '' [optional] Extra bowtie2 mapping options
	
	Paths can be absolute or relative.
	"
	}

# Defaults
CPUS=1
TMP=tmp
PAIREDEND=0

nopts=$#
for ((i=1 ; i <= nopts ; i++)); do
	case "$1" in
		-b | --bam)  # sorted BAM output file
			BAM="$2"
			echo "    output sorted and indexed BAM file: $2"
			shift 2
			;;
		-f | --fastq)  # single fastq input file
			IN_FASTQ="$2"
			echo "    input FastQ file: $2"
			shift 2
			;;
		-n | --numcpus)
			CPUS="$2"
			echo "    threads: $2"
			shift 2
			;;
		-pe | --pairedend)
			PAIREDEND="$2"
			echo "    paired-end mode: $2"
			shift 2
			;;
		-r | --ref)
			REFERENCE="$2"
			echo "    reference genome: $2"
			shift 2
			;;
		-x | --xopts)
			XOPTS="$2"
			echo "    extra options passed to bowtie2 for mapping: $2"
			shift 2
			;;
		-t | --temp)
			TMP="$2"
			echo "    temporary directory: $2"
			shift 2
			;;
		-v | --verbose)
			VERBOSE=1
			echo "    verbose mode invoked"
			shift
			;;
		-h | --help)
			usage
			exit 1
			;;
		\?)
			echo "    ERROR: $2 is not a valid argument"
			usage
			exit 1
			;;
	esac
done


# Input file requirements and dependency check
[[ -z "$IN_FASTQ" || -z "$REFERENCE" ]] && { usage; exit 1; }
[[ ! -s "$REFERENCE.4.bt2" || ! -s "$REFERENCE.rev.2.bt2" ]] && { echo 'ERROR: bowtie2 index file(s) not found'; exit 1; }
command -v bowtie2 >/dev/null 2>&1 || { echo 'ERROR: bowtie2 not found' >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo 'ERROR: samtools not found' >&2; exit 1; }


# Create tmp dir (if absent)
if [ ! -d "$TMP" ]; then
	mkdir -p "$TMP"
	echo '    Created temporary directory...'
else
	echo "WARNING: Temporary directory $TMP already exists. Files present in the specified tmp dir might disrupt analysis"
fi


# Deal with PE or SE reads
if [[ $PAIREDEND -eq 0 ]]; then
	# reassign var; script outputs 1 for PE and 0 for SE
	PAIREDEND=$((`run_assembly_isFastqPE.pl "$IN_FASTQ" --checkfirst 50` + 1))
	if [[ $PAIREDEND -eq 2 ]]; then
		echo "auto-detected $IN_FASTQ as paired-end..."
	elif [[ $PAIREDEND -eq 1 ]]; then
		echo "auto-detected $IN_FASTQ as single-end..."
	fi
fi

# Clip FASTQ[.GZ] and FQ[.GZ] extension (GNU sed required; not POSIX)
OUT_PREFIX=$(basename "$IN_FASTQ" | sed -r 's/\.(fastq|fq|fastq\.gz|fq\.gz)$//1')
if [[ $PAIREDEND -eq 2 ]]; then
	# bowtie2 requires separate files for PE reads
	run_assembly_shuffleReads.pl -d "$IN_FASTQ" 1>"$TMP"/"$OUT_PREFIX".1.fastq 2>"$TMP"/"$OUT_PREFIX".2.fastq
	BT2_READS="-1 $TMP/$OUT_PREFIX.1.fastq -2 $TMP/$OUT_PREFIX.2.fastq"
	echo 'deinterleaved PE file...'
elif [[ $PAIREDEND -eq 1 ]]; then
	# bowtie2 cannot accept gunzipped files
	if [ ${IN_FASTQ: -3} == '.gz' ]; then
		BT2_READS="-U <(zcat $IN_FASTQ)"
	else
		BT2_READS="-U $IN_FASTQ"
	fi
	echo 'using SE reads...'
else
	echo 'Usage: --pairedend <0|1|2> For auto, single-end, or paired-end respectively. Default: auto (0).'
	exit 1
fi


# Map with bowtie2
[[ $PAIREDEND -eq 2 ]] && BT2_READS="--no-discordant --no-mixed -X 1000 $BT2_READS"
COMMAND="$(echo bowtie2 -q -t -p $CPUS -k 10  \
	--sensitive-local -x $REFERENCE $BT2_READS \
	-S $TMP/$OUT_PREFIX.sam | sed 's/ \{1,\}/ /g')"
eval "$COMMAND"


# Post-processing of mapped sequences
if [ -s "$TMP"/"$OUT_PREFIX".sam ]; then
	# Convert output to binary, sort according to reference position/coordinate, and index
	samtools view -@ "$CPUS" -bSh "$TMP"/"$OUT_PREFIX".sam -o "$TMP"/"$OUT_PREFIX".unsorted.bam -T "$REFERENCE"
	samtools sort -@ "$CPUS" "$TMP"/"$OUT_PREFIX".unsorted.bam -O bam -o "$TMP"/"$OUT_PREFIX".sorted.bam -T "$BAM"
	mv -v "$TMP"/"$OUT_PREFIX".sorted.bam "$BAM"
	samtools index -b "$BAM"

	# Cleanup
	rm -fv "$TMP"/"$OUT_PREFIX".sam "$TMP"/"$OUT_PREFIX".unsorted.bam "$TMP"/"$OUT_PREFIX".[12].fastq
else
	echo ERROR: bowtie2 did not properly output the file "$TMP"/"$OUT_PREFIX".sam
	exit 1
fi

