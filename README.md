lyve-SET
========

LYVE version of the Snp Extraction Tool (SET), a method of using hqSNPs to create a phylogeny.  NML, part of PHAC, has the original version of SET.  However, I have been updating it since inception and so the results might differ slightly.

SET is meant to be run on a cluster but will just run system calls if qsub is not present. 

Requirements
------------
* Perl, multithreaded
* RAxML
* PhyML
* FreeBayes
* Smalt
* CG-Pipeline
* Samtools
* Schedule::SGELK (installed at the same time if you download with git)

Installation
------------
    $ git clone --recursive https://github.com/lskatz/lyve-SET.git

Usage
-----
    Usage: launch_set.pl -ref reference.fasta [-b bam/ -v vcf/ -t tmp/ -reads reads/ -m msa/]
      Where parameters with a / are directories
      -r where fastq and fastq.gz files are located
      -b where to put bams
      -v where to put vcfs
      --msadir multiple sequence alignment and tree files (final output)
      -numcpus number of cpus
      -numnodes maximum number of nodes
      -w working directory where qsub commands can be stored. Default: CWD
      -a allowed flanking distance in bp. Nucleotides this close together cannot be considered as high-quality.
      --nomsa to not make a multiple sequence alignment
      --notrees to not make phylogenies
      -q '-q long.q' extra options to pass to qsub. This is not sanitized.
      --noclean to not clean reads before mapping (faster, but you need to have clean reads to start with)
