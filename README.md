lyve-SET
========

LYVE version of the Snp Extraction Tool (SET), a method of using hqSNPs to create a phylogeny.  NML, part of PHAC, has the original version of SET (https://github.com/apetkau).  However, I have been updating it since inception and so the results might differ slightly.

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
* wgsim (part of samtools) if you are simulating reads
* Schedule::SGELK (installed at the same time if you download with git)
* CallSam (installed at the same time if you download with git)

Installation
------------

To install, just run the make file with `make install`.  Run `make help` for all options.

Usage
-----
    Usage: launch_set.pl -ref reference.fasta [-b bam/ -v vcf/ -t tmp/ -reads reads/ -m msa/ -asm asm/]
    Where parameters with a / are directories
    -reads    readsdir/       where fastq and fastq.gz files are located
    -bam      bamdir/         where to put bams
    -vcf      vcfdir/         where to put vcfs
    --tmpdir  tmpdir/         tmp/ Where to put temporary files
    --msadir  msadir/         multiple sequence alignment and tree files (final output)
    -asm      asmdir/         directory of assemblies. Copy or symlink the reference genome assembly to use it if it is not already in the raw reads directory
    -all      0               allowed flanking distance in bp. Nucleotides this close together cannot be considered as high-quality.
      NOTE: Set -all to 'auto' to let SET determine this distance using snpDistribution.pl

    SKIP CERTAIN STEPS
    --noclean to not clean reads before mapping (faster, but you need to have clean reads to start with; removes the requirement for CG-Pipeline)
    --nomsa to not make a multiple sequence alignment
    --notrees to not make phylogenies
    MODULES
    --mapper       smalt             Which mapper? Choices: smalt, snap
    --snpcaller    freebayes         Which SNP caller? Choices: freebayes, callsam
    --msa-creation lyve-set          Which method of making the multiple sequence alignment? lyve-set, lyve-set-lowmem (unvalidated)
    SCHEDULER AND MULTITHREADING OPTIONS
    --queue     all.q         The default queue to use.
    --qsubxopts '-N lyve-set' extra options to pass to qsub. This is not sanitized; internal options might overwrite yours.
    --numnodes  20  maximum number of nodes
    --numcpus   1  number of cpus
    -w dir/     working directory where qsub commands can be stored. Default: CWD/.SGELK/


Examples
------
If you have version 0.8 or earlier, use the following sytax:

    # Set up the directory with your reads and reference genome
    $ set_create.pl setTest  # directory structure created
    $ tree setTest           # View the new directory structure
    setTest
    ├── asm
    ├── bam
    ├── msa
    ├── reads
    ├── reference
    ├── tmp
    └── vcf
        └── unfiltered

    $ cd setTest/            
    # copy over all your fastq files to the reads directory
    $ cp -v ../path/to/fastq.gz/dir/*.fastq.gz reads/     
    # OR symlink your fastq files
    $ cd reads
    $ ln -sv ../../path/to/fastq.gz/dir/*.fastq.gz .
    $ cd ..
    # copy over all assembled genomes that you want to include in the tree in the asm directory
    $ cp -v ../path/to/fasta/dir/*.fasta asm/                 
    # copy over your single reference genome (only one permitted)
    $ cp -v ../path/to/fasta/dir/reference.fasta reference/  

If you have a newer version >0.8 then you have the script set_manage.pl and you should use the following syntax:
    
    $ set_manage.pl -c setTest
    $ set_manage.pl setTest --add-reads file1.fastq.gz
    $ set_manage.pl setTest --add-reads file2.fastq.gz
    $ set_manage.pl setTest --add-reads file3.fastq.gz
    $ set_manage.pl setTest --add-reads file4.fastq.gz
    $ set_manage.pl setTest --add-reads file5.fastq.gz
    $ set_manage.pl setTest --add-assembly file1.fasta
    $ set_manage.pl setTest --add-assembly file2.fasta
    $ set_manage.pl setTest --change-reference file3.fasta
    $ cd  setTest # get into the directory before running launch_set.pl

NOTE: no underscores or dashes allowed in the reference genome fasta file
    
Run Lyve-SET with as few options as possible

    $ cd setProj && launch_set.pl -ref reference/reference.fasta  # < v0.8.1
    $ launch_set.pl setProj                                       # >= v0.8.1

More complex

    # < v0.8.1
    $ cd setProj && launch_set.pl -ref reference/reference.fasta  --queue all.q --numnodes 20 --numcpus 16 --noclean --notrees
    # >= v0.8.1
    $ launch_set.pl setProj --queue all.q --numnodes 20 --numcpus 16 --noclean --notrees
    # Change some module steps, >= v0.7
    $ launch_set setProj --snpcaller callsam --msa-creation lyve-set-lowmem
    
If you specified notrees, then you can edit the multiple sequence alignment before analyzing it. See the next section on examples on how/why you would edit the alignment.

    $ cd msa
    $ gedit out.aln.fas  # alter the deflines or whatever you want before moving on
    # => out.aln.fas is here
    $ set_process_msa.pl --auto --numcpus 12
    # Optionally, qsub this script instead because it could be cpu-intensive
    $ qsub -pe smp 12 -cwd -V -o trees.log -j y -N msaLyveSET -S $(which perl) $(which set_process_msa.pl) --auto --numcpus 12

Why would you want to edit the out.aln.fas file?  Or what kinds of things can you observe here before making a tree?
    
    # Alter the identifiers of your genomes, so that they look nice in the phylogeny(ies)
    $ sed -i.bak 's/\.fastq\.gz.*//' out.aln.fas

    # After the taxon names are edited nicely,
    # make a new file from which you can extract your taxa of choice (out2.aln.fas).
    $ cp -v out.aln.fas out2.aln.fas

    # => All extensions are removed in taxon names; a backup of the file was named out.aln.fas.bak
    # find the genomes with the most number of Ns (ie masked SNP calls)
    $ perl -lane 'chomp; if(/^>/){s/>//;$id=$_;}else{$num=(s/(N)/$1/gi); print "$id\t$num";}' < out2.aln.fas|sort -k2,2n|column -t
    # => consider removing any genome with too many masked bases

    # Run SET on your new set of genomes out2.aln.fas.
    $ set_process_msa.pl out2.aln.fas --auto --numcpus 12
    
    # Create a new subset as needed, but you should always read from your master record,
    # which is the original out.aln.fas.
    $ cp -v out.aln.fas out3.aln.fas
    $ ...
    $ set_process_msa.pl out3.aln.fas --auto --numcpus 12

Citing lyve-SET
-----
To cite lyve-SET, please reference this site and cite the Haiti Anniversary paper. Lyve-SET also makes use of the tools shown above in the prerequisites.  If you feel like your study relied heavily on any of those tools, please don't forget to cite them!
    
    https://github.com/lskatz/lyve-SET
    Katz LS, Petkau A, Beaulaurier J, Tyler S, Antonova ES, Turnsek MA, Guo Y, Wang S, Paxinos EE, Orata F, Gladney LM, Stroika S, Folster JP, Rowe L, Freeman MM, Knox N, Frace M, Boncy J, Graham M, Hammer BK, Boucher Y, Bashir A, Hanage WP, Van Domselaar G, Tarr CL. 
    2013. Evolutionary dynamics of Vibrio cholerae O1 following a single-source introduction to Haiti. 
    MBio 4(4):e00398-13. doi:10.1128/mBio.00398-13.
