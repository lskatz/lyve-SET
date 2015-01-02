Examples
========

Run a test dataset
------------------

The script `set_test.pl` will run an actual test on a given dataset. It uses `set_downloadTestData.pl` to get any bacterial genomes and then runs `launch_set.pl`.  However, the lambda dataset is small enough to fit on GIT and does not need to be downloaded.

    Runs a test dataset with Lyve-SET
    Usage: set_test.pl dataset [project]
    dataset names could be one of the following:
      escherichia_coli, lambda, listeria_monocytogenes, salmonella_enterica_agona
    NOTE: project will be the name of the dataset, if it is not given

    --numcpus 1  How many cpus you want to use
    --do-nothing To print the commands but do not run system calls

`$ set_test.pl lambda  # will run the entire lambda phage dataset and produce meaningful results in ./lambda/msa/`


Prepare the project directory
-----------------------------

The script `set_manage.pl` sets up the project directory and adds reads, and you should use the following syntax:
    
    $ set_manage.pl --create setTest

Depending on your knowledge of Linux, you might choose to set up the rest of the project using `set_manage.pl` or using symlinks.  This is the `set_manage.pl` way:

    $ set_manage.pl setTest --add-reads file1.fastq.gz
    $ set_manage.pl setTest --add-reads file2.fastq.gz
    $ set_manage.pl setTest --add-reads file3.fastq.gz
    $ set_manage.pl setTest --add-reads file4.fastq.gz
    $ set_manage.pl setTest --add-reads file5.fastq.gz
    $ set_manage.pl setTest --add-assembly file1.fasta
    $ set_manage.pl setTest --add-assembly file2.fasta
    $ set_manage.pl setTest --change-reference file3.fasta

This is the symlink way:

    $ cd setTest/reads
    $ ln -sv path/to/reads/*.fastq.gz .   # symlink reads
    $ cd ../asm
    $ ln -sv path/to/assemblies/*.fasta . # symlink the assemblies
    $ cd ../reference
    # copy the assembly in case you need to alter it (e.g., remove small contigs or edit deflines)
    $ cp -v path/to/assemblies/reference.fasta .

Run Lyve-SET with as few options as possible

    $ launch_set.pl setTest

More complex

    $ launch_set.pl setTest --queue all.q --numnodes 20 --numcpus 12 --noclean --notrees
    
If you specified notrees, then you can edit the multiple sequence alignment before analyzing it. See the next section on examples on how/why you would edit the alignment.

    $ cd setTest/msa
    $ gedit out.aln.fas  # alter the deflines or whatever you want before moving on
    # => out.aln.fas is here
    $ set_process_msa.pl --auto --numcpus 12
    # Optionally, qsub this script instead because it could be cpu-intensive
    $ qsub -pe smp 12 -cwd -V -o trees.log -j y -N msaLyveSET -S $(which perl) $(which set_process_msa.pl) --auto --numcpus 12

If you make a mistake and need to redo something:

    # remove all intermediate files
    $ rm setProj/bam/genome*
    $ rm setProj/vcf/genome* setProj/vcf/unfiltered/genome*
    # OR, remove a genome entirely
    $ set_manage.pl setProj --remove-reads genome.fastq.gz
    $ set_manage.pl setProj --remove-assembly genome.fasta
    # remove the last multiple sequence alignment files
    $ rm -r setProj/msa/*
    # or save the MSA results for another time
    $ mv setProj/msa setProj/msa.oldresults && mkdir setProj/msa
    # redo the analysis (all untouched bams and vcfs will not be redone)
    $ launch_set.pl setProj

Why would you want to edit the out.aln.fas file?  Or what kinds of things can you observe here before making a tree?
    
    # Alter the identifiers of your genomes, so that they look nice in the phylogeny(ies)
    $ sed -i.bak 's/\.fastq\.gz.*//' out.aln.fas
    # OR maybe to just remove the reference genome name:
    $ sed -i.bak 's/\-reference.*//' out.aln.fas

    # Be sure that the taxa names are unique still.
    # If there is any output, then you have duplicated names which need to be fixed.
    $ grep ">" out.aln.fas | sort | uniq -d # nothing should show up

    # After the taxon names are edited nicely,
    # back up the out.aln.fas file before any entries are removed
    $ cp -v out.aln.fas out.aln.fas.bak

    # => All extensions are removed in taxon names; a backup of the file was named out.aln.fas.bak
    # find the genomes with the most number of Ns (ie masked SNP calls)
    $ perl -lane 'chomp; if(/^>/){s/>//;$id=$_;}else{$num=(s/(N)/$1/gi); print "$id\t$num";}' < out.aln.fas|sort -k2,2n|column -t
    # => consider removing any genome with too many masked bases by manually editing the file.

    # Run SET on your new set of genomes out.aln.fas.
    $ set_process_msa.pl out.aln.fas --auto --numcpus 12
    
    # Create a new subset as needed, but you should always read from your master record and not the informative.aln.fas file
    # which is now out.aln.fas.bak
    $ cp -v out.aln.fas.bak out.aln.fas
    $ ... # more editing here
    $ set_process_msa.pl out.aln.fas --auto --numcpus 12


Specific ways to regenerate files
---------------------------------

Lyve-SET is very modular and so there are specific scripts to regenerate files.  Since you might edit the msa files, you might want to know how to recover in case you make a mistake.

Need to remake out.aln.fas?

    $ mvcfToAlignment.pl out.pooled.vcf.gz --bcfOutput bcfquery.out > out.aln.fas

Need to remake out.pooled.vcf.gz? Use bcftools.

    $ bcftools merge vcf/unfiltered/*.vcf.gz -O -z > msa/out.pooled.vcf.gz
    $ tabix -f msa/out.pooled.vcf.gz # Always index your compressed vcf files

Need to remake all the other msa files after you recreated out.aln.fas?
  
    $ set_processMsa.pl --auto out.aln.fas --numcpus 12 --force

