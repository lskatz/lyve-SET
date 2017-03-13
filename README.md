Lyve-SET
========

LYVE version of the Snp Extraction Tool (SET), a method of using hqSNPs to create a phylogeny.  Lyve-SET is meant to be run on a cluster but will just run system calls if qsub is not present.

NML, part of PHAC, has the original version of SET (https://github.com/apetkau).  However, I have been updating it since inception and so the results might differ slightly.

Installation
------------
* `make install`
* `make help` - for other `make` options
* See [INSTALL.md](docs/INSTALL.md) for more information including prerequisite software

For the impatient
-----------------
Here is a way to just try out the test dataset. It runs the lambda dataset (1st param) and places the results into the lambda directory (2nd param)

    set_test.pl --numcpus 8 lambda lambda
    # try out a larger dataset
    set_test.pl --numcpus 8 listeria_monocytogenes listeria_monocytogenes 

Make Lyve-SET go quickly with `--fast`!  This option is shorthand for several other options that save on computational time. See `launch_set.pl` usage below for more details.

    set_test.pl --numcpus 8 --fast listeria_monocytogenes listeria_monocytogenes 

Usage
-----
To see the help for any script, run it without options or with `--help`.  For example, `set_test.pl -h`.  The following is the help for the main script, `launch_set.pl`:

    Usage: launch_set.pl [project] [-ref reference.fasta|reference.gbk]

    If project is not given, then it is assumed to be the current working directory.

    -ref               ref.fasta      The reference genome assembly. If it is
                                      a genbank or embl file, then it will be
                                      converted to reference.gbk.fasta and will
                                      be used for SNP annotation. If a fasta
                                      is given, then no SNP annotation will
                                      happen.
                                      Default: project/reference/reference.fasta

    COMMON OPTIONS
    --allowedFlanking  0              allowed flanking distance in bp.
                                      Nucleotides this close together cannot be
                                      considered as high-quality.
    --min_alt_frac     0.75           The percent consensus that needs
                                      to be reached before a SNP is called.
                                      Otherwise, 'N'
    --min_coverage     10             Minimum coverage needed before a
                                      SNP is called. Otherwise, 'N'
    --presets          ""             See presets.conf for more information
    --numcpus          1              number of cpus

    PERFORM CERTAIN STEPS
    --mask-phages                  Search for and mask phages in the reference genome
    --mask-cliffs                  Search for and mask 'Cliffs' in pileups

    SKIP CERTAIN STEPS
    --nomatrix                     Do not create an hqSNP matrix
    --nomsa                        Do not make a multiple sequence alignment
    --notrees                      Do not make phylogenies
    --singleend                    Treat everything like single-end. Useful
                                   for when you think there is a single-
                                   end/paired-end bias.
    OTHER SHORTCUTS
    --fast                         Shorthand for --downsample --mapper snap --nomask-phages
                                                 --nomask-cliffs --sample-sites
    --downsample                   Downsample all reads to 50x. Approximated according
                                   to the ref genome assembly
    --sample-sites                 Randomly choose a genome and find SNPs in a quick
                                   and dirty way. Then on the SNP-calling stage,
                                   only interrogate those sites for SNPs for each
                                   genome (including the randomly-sampled genome).

    MODULES
    --read_cleaner none            Which read cleaner? Choices: none, CGP, BayesHammer
    --mapper       smalt           Which mapper? Choices: smalt, snap
    --snpcaller    varscan         Which SNP caller? Choices: varscan, vcftools

    SCHEDULER AND MULTITHREADING OPTIONS
    --queue        all.q           default queue to use
    --numnodes     50              maximum number of nodes
    --qsubxopts    '-N lyve-set'   Extra options to pass to qsub. This is not
                                   sanitized; internal options might overwrite yours.
    --noqsub                       Do not use the scheduler, even if it exists

    LOCATIONS OF FILE DIRECTORIES
    --readsdir  readsdir/          where fastq and fastq.gz files are located
    --bamdir    bamdir/            where to put bams
    --vcfdir    vcfdir/            where to put vcfs
    --tmpdir    tmpdir/            tmp/ Where to put temporary files
    --msadir    msadir/            multiple sequence alignment and tree files (final output)
    --logdir    logdir/            Where to put log files. Qsub commands are also stored here.
    --asmdir    asmdir/            directory of assemblies. Copy or symlink the reference genome assembly
                                   to use it if it is not already in the raw reads directory

Run a test dataset
------------------

See: [examples.md](docs/EXAMPLES.md) for more details.  
Also see: [testdata.md](docs/TESTDATA.md) for more details on making your own test data set.

The script `set_test.pl` will run an actual test on a given dataset

    Runs a test dataset with Lyve-SET
    Usage: set_test.pl dataset [project]
    dataset names could be one of the following:
      escherichia_coli, lambda, listeria_monocytogenes, salmonella_enterica_agona
    NOTE: project will be the name of the dataset, if it is not given

    --numcpus 1  How many cpus you want to use
    --do-nothing To print the commands but do not run system calls

    # will run the entire lambda phage dataset and produce meaningful results in ./lambda/msa/
    $ set_test.pl lambda lambda


Examples
------

See: [examples.md](docs/EXAMPLES.md) for more details.

The script `set_manage.pl` sets up the project directory and adds reads, and you should use the following syntax. Note that paired end reads should be in interleaved format. Scripts that interleave reads include `run_assembly_shuffleReads.pl` in the CG-Pipeline package (included with `make install`) and also `shuffleSequences_fastq.pl` in the Velvet package.
    
    # Shuffle your reads if they are not already.
    $ mkdir interleaved
    $ for f in *_R1.fastq.gz; do
    >   b=$(basename $i _R1.fastq.gz)  # getting the basename of the file
    >   r=${b}_R2.fastq.gz             # Reverse reads filename
    >   run_assembly_shuffleReads.pl $f $r | gzip -c > interleaved/$b.fastq.gz
    > done;
    # Create the project directory `setTest`
    $ set_manage.pl --create setTest
    # Add reads
    $ for i in interleaved/*.fastq.gz; do
    >   set_manage.pl setTest --add-reads $i
    > done;
    # Add assemblies (optional)
    $ set_manage.pl setTest --add-assembly file1.fasta
    $ set_manage.pl setTest --add-assembly file2.fasta
    # Specify your reference genome
    $ set_manage.pl setTest --change-reference file3.fasta

    
Run Lyve-SET with as few options as possible

    $ launch_set.pl setProj

More complex

    $ launch_set.pl setProj --queue all.q --numnodes 20 --numcpus 16 --notrees
    
Output files
------------
Most output files that you will want to see are under project/msa.  However for more details please see [docs/output.md](docs/OUTPUT.md).

Getting help
------------
* Check out the [FAQ](docs/FAQ.md) first to see if your question has already been asked
* Join the [Google Group](https://groups.google.com/forum/#!forum/lyve-set) at https://groups.google.com/forum/#!forum/lyve-set 
* [Tips and tricks](docs/TIPS.md)

Citing lyve-SET
-----
    
    https://github.com/lskatz/lyve-SET
    Katz LS, Griswold T, Williams-Newkirk AJ, Wagner D, Petkau A, et al. (2017) A Comparative Analysis of the Lyve-SET Phylogenomics Pipeline for Genomic Epidemiology of Foodborne Pathogens. Frontiers in Microbiology 8.
    
