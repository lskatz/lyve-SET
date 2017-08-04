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
  * You can use the [lskScript](https://github.com/lskatz/lskScripts) tool `validateFastq.pl` to find common mistakes.  For example: `zcat reads/genome.fastq.gz | validateFastq.pl --pe --min-length 1 --verbose`.
  * Is the interleaved file line count the number of lines of both R1 and R2? Example command: `zcat R1.fastq.gz | wc -l`.
  * Simply reshuffle: `shuffleSplitReads.pl -o shuffled --numcpus 1 split/*.fastq.gz`

## BCFtools merge

### cannot read fastq header; error with bcftools merge

This is likely due to [Fastq integrity](#fastq-integrity)

## The tree was not created

### Too few taxa

RAxML requires at least four genomes to create a tree.  If you have few genomes in your analysis, perhaps the `out.pairwise.tsv` or `out.pairwiseMatrix.tsv` file would make more sense in terms of deliverable results.

### Nonunique names

* RAxML only considers the **first 30 characters** of your genome name.  In the offhand chance that you have two genomes with filenames that vary on the 31st or later characters, then consider renaming your files appropriately.
* Did you include **R1 and R2 and/or also the shuffled reads**?  Lyve-SET requires shuffled reads and so if you haven't done this, please follow the steps under the section [Interleaved reads](#interleaved-reads).  You should not be using both R1 and R2 in your input directly.  If they are already shuffled, delete all instances of R1 and R2 files in the project directory.  However, *do not* delete the shuffled file.  *Do not* delete your original data (just a general tip).  Delete R1/R2 files and their derivative files in the reads subdirectory, the tmp subdirectory, the bam subdirectory, and the vcf subdirectory.  Delete all files in the msa directory, because they will have to be recreated.  Rerun Lyve-SET.  This is a huge issue because
  1. R1 and R2 SNPs will not be as accurate as the paired-end shuffled entry, because their read mapping is not as accurate
  2. These less accurate entries will enter noise into the tree inferrence.

### Not enough variable sites

This isn't a huge problem but it's not immediately obvious that this is the case always.  If all of your genomes have the exact same SNP profile, then the multiple sequence alignment will have zero variable sites, and RAxML will fail.  This is a result in itself because it says that all your genomes have the same genomic profile.

How do you detect this?  You can view `out.snpmatrix.tsv` to see if there are many sites without majority ambiguous allele calls.  You can see if `out.snpmatrix.tsv` has many sites.  Therefore you can do a simple calculation to show that you have many high-quality sites but few or no variable sites.

How do you avoid this situation?  If you have a good outgroup, then there will be variable sites, and it won't matter if your clade of interest has no variation between them.  RAxML can build a tree when there are a few guaranteed variable sites introduced by the outgroup.  This outgroup by definition is a genome that is phylogenetically related but is not the same profile as your clade of interest.
