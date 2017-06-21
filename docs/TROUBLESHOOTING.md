# Troubleshooting

Things do not always go smoothly.  Here is a list of things to try out if things go wrong.

## Interleaved reads

Why are both R1 and R2 showing up in my results?

Lyve-SET requires _interleaved_ reads (or commonly called "shuffled").  This means that there is one fastq file per genome and that their format is such that there is one forward read, followed by its reverse reads, followed by the next forward read, etc.
To interleave one set of reads, you can use the script `run_assembly_shuffleReads.pl`.  To shuffle many pairs of reads, you can use `shuffleSplitReads.pl`.  Both have usage statements if you run them with `--help`.

## Smalt says that there is an invalid fastq/fasta format

This can be the result of a few different things.  

### Reference assembly integrity

First, check to see if your reference genome is intact.

* Manually inspect the fasta file with `less`, e.g., `less reference/reference.fasta`.
* Another way to manually inspect the fasta file is with `grep`, e.g., `grep -A 1 ">" reference/reference.fasta`.
* Get some summary metrics with `run_assembly_metrics.pl`, e.g., `run_assembly_metrics.pl reference/reference.fasta | column -t`.

### Fastq integrity

You can also see if your fastq files are intact which is a little more tricky.

* Are your files in the interleaved format? See [Interleaved reads](#interleaved-reads).
* Are your files gzipped?  Use the Linux command `file` to see that it says 'gzip' in the description of your files.  For example, `file reads/*.fastq.gz`.  It is not necessary to have gzipped files, but it _is_ necessary for the extension to match the type. For example you do not want to have compressed files that end in `.fastq.gz` or uncompressed files that end in `.fastq`.
* Are the actual files intact?  
** You can use the [lskScript](https://github.com/lskatz/lskScripts) tool `validateFastq.pl` to find common mistakes.  For example: `zcat reads/genome.fastq.gz | validateFastq.pl --pe --min-length 1 --verbose`.
** Is the interleaved file line count the number of lines of both R1 and R2? Example command: `zcat R1.fastq.gz | wc -l`.
** Simply reshuffle: `shuffleSplitReads.pl -o shuffled --numcpus 1 split/*.fastq.gz`

## BCFtools merge

### cannot read fastq header; error with bcftools merge

This is likely due to [# Fastq integrity]

