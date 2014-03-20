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
* Schedule::SGELK (installed at the same time if you download with git)

Installation
------------
    
    $ mkdir ~/bin
    $ cd ~/bin
    $ git clone --recursive https://github.com/lskatz/lyve-SET.git
    $ export PATH=$PATH:~/bin/lyve-SET # you might also put this line into your .bash_profile or other login script

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
Lyve-SET is modular and so the individual scripts can be run too.  For example, you can run launch\_smalt.pl or launch\_snap.pl to run mapping alone; however, indexing the reference fasta takes place in launch_set.pl.  To get usage help on any of these scripts, run the script with no options.

Citing lyve-SET
-----
To cite lyve-SET, please reference this site and cite the Haiti Anniversary paper. Lyve-SET also makes use of the tools shown above in the prerequisites.  If you feel like your study relied heavily on any of those tools, please don't forget to cite them!
    
    https://github.com/lskatz/lyve-SET
    Katz LS, Petkau A, Beaulaurier J, Tyler S, Antonova ES, Turnsek MA, Guo Y, Wang S, Paxinos EE, Orata F, Gladney LM, Stroika S, Folster JP, Rowe L, Freeman MM, Knox N, Frace M, Boncy J, Graham M, Hammer BK, Boucher Y, Bashir A, Hanage WP, Van Domselaar G, Tarr CL. 2013. Evolutionary dynamics of Vibrio cholerae O1 following a single-source introduction to Haiti. mBio 4(4):e00398-13. doi:10.1128/mBio.00398-13.
