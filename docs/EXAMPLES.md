Examples
========

Run a test dataset
------------------

See: [examples.md](docs/EXAMPLES.md) for more details.

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

NOTE: no underscores or dashes allowed in the reference genome fasta file headers
    
Run Lyve-SET with as few options as possible

    $ launch_set.pl setProj

More complex

    $ launch_set.pl setProj --queue all.q --numnodes 20 --numcpus 16 --noclean --notrees
    
If you specified notrees, then you can edit the multiple sequence alignment before analyzing it. See the next section on examples on how/why you would edit the alignment.

    $ cd setProj/msa
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

