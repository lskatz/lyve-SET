lyve-SET
========

LYVE version of the Snp Extraction Tool (SET), a method of using hqSNPs to create a phylogeny.  Lyve-SET is meant to be run on a cluster but will just run system calls if qsub is not present.

NML, part of PHAC, has the original version of SET (https://github.com/apetkau).  However, I have been updating it since inception and so the results might differ slightly.

Requirements
------------
* **Must-have** and not installed with `make install`
  * **Perl, multithreaded**
  * **RAxML**
  * **Smalt** (in latest release but not cutting edge version)
* **Requirements installed with `make install`**
  * CG-Pipeline
  * Schedule::SGELK
  * VCF Tools and Vcf.pm
  * Samtools
    * wgsim if you are simulating reads
    * bgzip
    * tabix
  * PHAST
  * Varscan
* Optional "requirements"
  * PhyML
  * FreeBayes
  * SNAP

Installation
------------
* `make install`
* `make check` - check and see if you have all the prerequisites
* `make test` - run a test phage dataset provided by CFSAN
* `make help` - for other `make` options

For the impatient
-----------------
Here is a way to just try out the test dataset.
    set_test.pl lambda --numcpus 8 # or however many cpus you want
    set_test.pl listeria_monocytogenes --numcpus 8 # or another dataset

Usage
-----
    Usage: launch_set.pl [project] [-ref reference.fasta]
    If project is not given, then it is assumed to be the current working directory.
    If reference is not given, then it is assumed to be proj/reference/reference.fasta
    Where parameters with a / are directories
    -ref      proj/reference/reference.fasta   The reference genome assembly
    -reads    readsdir/       where fastq and fastq.gz files are located
    -bam      bamdir/         where to put bams
    -vcf      vcfdir/         where to put vcfs
    --tmpdir  tmpdir/         tmp/ Where to put temporary files
    --msadir  msadir/         multiple sequence alignment and tree files (final output)
    --logdir  logdir/         Where to put log files. Qsub commands are also stored here.
    -asm      asmdir/         directory of assemblies. Copy or symlink the reference genome assembly to use it if it is not already in the raw reads directory

    SNP MATRIX OPTIONS
    --allowedFlanking  0               allowed flanking distance in bp. Nucleotides this close together cannot be considered as high-quality.  Set to -1 to let SET determine this distance using snpDistribution.pl
    --min_alt_frac     0.75  The percent consensus that needs to be reached before a SNP is called. Otherwise, 'N'
    --min_coverage     10  Minimum coverage needed before a SNP is called. Otherwise, 'N'

    SKIP CERTAIN STEPS
    --noclean to not clean reads before mapping (faster, but you need to have clean reads to start with; removes the requirement for CG-Pipeline)
    --nomatrix to not create an hqSNP matrix
    --nomsa to not make a multiple sequence alignment
    --notrees to not make phylogenies
    MODULES
    --mapper       smalt             Which mapper? Choices: smalt, snap
    --snpcaller    freebayes         Which SNP caller? Choices: freebayes, varscan
    SCHEDULER AND MULTITHREADING OPTIONS
    --queue     all.q         The default queue to use.
    --qsubxopts '-N lyve-set' extra options to pass to qsub. This is not sanitized; internal options might overwrite yours.
    --numnodes  20  maximum number of nodes
    --numcpus   1  number of cpus

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

`$ set_test.pl lambda  # will run the entire lambda phage dataset and produce meaningful results in ./lambda/msa/`


Examples
------

See: [examples.md](docs/EXAMPLES.md) for more details.

The script `set_manage.pl` sets up the project directory and adds reads, and you should use the following syntax:
    
    $ set_manage.pl --create setTest
    $ set_manage.pl setTest --add-reads file1.fastq.gz
    $ set_manage.pl setTest --add-reads file2.fastq.gz
    $ set_manage.pl setTest --add-reads file3.fastq.gz
    $ set_manage.pl setTest --add-reads file4.fastq.gz
    $ set_manage.pl setTest --add-reads file5.fastq.gz
    $ set_manage.pl setTest --add-assembly file1.fasta
    $ set_manage.pl setTest --add-assembly file2.fasta
    $ set_manage.pl setTest --change-reference file3.fasta

NOTE: paired end reads should be in interleaved format. Scripts that interleaved include run_assembly_shuffleReads.pl in the CG-Pipleline package and also shuffleSequences_fastq.pl in the Velvet package.

Run Lyve-SET with as few options as possible

    $ launch_set.pl setProj

More complex

    $ launch_set.pl setProj --queue all.q --numnodes 20 --numcpus 16 --noclean --notrees
    
Output files
------------
Most output files that you will want to see are under project/msa.  However for more details please see [docs/output.md](docs/OUTPUT.md).

Citing lyve-SET
-----
To cite lyve-SET, please reference this site and cite the Haiti Anniversary paper. Lyve-SET also makes use of the tools shown above in the prerequisites.  If you feel like your study relied heavily on any of those tools, please don't forget to cite them!
    
    https://github.com/lskatz/lyve-SET
    Katz LS, Petkau A, Beaulaurier J, Tyler S, Antonova ES, Turnsek MA, Guo Y, Wang S, Paxinos EE, Orata F, Gladney LM, Stroika S, Folster JP, Rowe L, Freeman MM, Knox N, Frace M, Boncy J, Graham M, Hammer BK, Boucher Y, Bashir A, Hanage WP, Van Domselaar G, Tarr CL. 
    2013. Evolutionary dynamics of Vibrio cholerae O1 following a single-source introduction to Haiti. 
    MBio 4(4):e00398-13. doi:10.1128/mBio.00398-13.
