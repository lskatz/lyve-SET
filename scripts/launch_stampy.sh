#!/bin/bash

# Usage: $0 -f file.fastq -b file.bam -t tmp/ -r reference.fasta
  # -t tmp to set the temporary directory as 'tmp'
  # --numcpus 1 number of cpus to use
  # -s '' Extra smalt map options (not validated). Default: $$settings{smaltxopts} 
  # --pairedend <0|1|2> For 'auto', single-end, or paired-end respectively. Default: auto (0).
  # --minPercentIdentity 95  The percent identity between a read and its match before it can be mapped

SCRIPT=$(basename $0);
function usage() { 
	echo "
	Usage: $SCRIPT -f file.fastq[.gz] -b file.bam -t tmp/ -r reference.fasta [-s '']
		
		-t [default: ./tmp] Path to temporary directory
		-c [default: 1] Number of threads to use
		-pe <0|1|2> [default: 0] For 'auto', single-end, or paired-end respectively
		
		-s '' [optional] Extra stampy mapping options
	
	Paths can be absolute or relative.
	"
	}

# Defaults
CPUS=1
TMP=tmp
PAIREDEND=0

nopts=$#
for ((i=1 ; i <= nopts ; i++)); do
  echo "$i $1 $2";
	case "$1" in
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
		-s | --xopts)
			XOPTS="$2"
			echo "    extra options passed to stampy for mapping: $2"
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

if [ "$IN_FASTQ" == "" ]; then
  usage
  exit 1
fi

# Create tmp dir (if absent)
if [ ! -d "$TMP" ]; then
	mkdir -p "$TMP"
	echo '    Created temporary directory...'
else
	echo "WARNING: temporary directory $TMP already exists. Files present in the specified tmp dir might disrupt analysis"
fi

OUT_PREFIX=$(basename "$IN_FASTQ" .fastq)

# Deal with PE or SE reads
if [[ $PAIREDEND -eq 0 ]]; then
	# reassign var; script outputs 1 for PE and 0 for SE
	PAIREDEND=$((`run_assembly_isFastqPE.pl "$IN_FASTQ" --checkfirst 50` + 1))
	continue
elif [[ $PAIREDEND -eq 2 ]]; then
	# stampy requires separate files for PE reads
	run_assembly_shuffleReads.pl -d "$IN_FASTQ" 1>"$OUT_PREFIX".1.fastq 2>"$OUT_PREFIX".2.fastq
	STAMPY_READS=$(-M "$OUT_PREFIX".1.fastq,"$OUT_PREFIX".2.fastq)
	echo 'deinterleaved PE file...'
elif [[ $PAIREDEND -eq 1 ]]; then
	STAMPY_READS=$(-M "$IN_FASTQ")
	echo 'using SE reads...'
else
	echo 'Usage: --pairedend <0|1|2> For auto, single-end, or paired-end respectively. Default: auto (0).'
	exit 1
fi

# Create reference index and hash table
stampy.py -G "$REFERENCE" "$TMP"/"$REFERENCE" 
stampy.py -g "$TMP"/"$REFERENCE".stidx -H "$TMP"/"$REFERENCE"

# Map with stampy
stampy.py -g "$TMP"/"$REFERENCE".stidx -h "$TMP"/"$REFERENCE".sthash \
--minposterior=20 --sensitive -t "$CPUS" \
--inputformat=fastq -M "$STAMPY_READS" -f sam \
--insertsize=300 --insertsd=250 -o "$TMP"/"OUT_PREFIX".sam "$XOPTS"

# Convert output to binary, sort according to reference position/coordinate, and index
samtools view -@ "$CPUS" -bSh "$TMP"/"$OUT_PREFIX".sam -o "$TMP"/"$OUT_PREFIX"_unsorted.bam
samtools sort -@ "$CPUS" "$TMP"/"$OUT_PREFIX"_unsorted.bam -O bam -o "$TMP"/"$OUT_PREFIX".sorted.bam
samtools index -b "$TMP"/"$OUT_PREFIX".sorted.bam


# Cleanup
rm -fv "$TMP"/"REFERENCE".stidx "$TMP"/"$REFERENCE".stidx
rm -fv "$TMP"/"OUT_PREFIX".sam "$TMP"/"$OUT_PREFIX"_unsorted.bam

