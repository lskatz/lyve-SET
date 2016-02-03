#Reference genomes

##Can I have multiple reference genomes?

No.  However, you can add multiple assemblies to the asm directory.  Assemblies don't have the same error profile as reads and so you might expect some skew in the ultimate phylogeny.

##How do I choose a reference genome?

The best reference genome for an outbreak is something in-clade.  You might even need to assemble the genome from your own reads before starting.  An outgroup is not as related to your clade by definition and so it is probably not the best genome to use.  The next best quality is a closed genome, or a genome with a high N50.

#My custom input files

##How do I include an SFF file?

There currently is no way to use an SFF file directly.  These usually come from 454 or Ion Torrent and therefore have a different error profile than Illumina.  Therefore, you might observe some bias if you mix these chemistries.  

* You can convert your SFF files to fastq using a variety of tools.  Some of these include `sffinfo` from the Newbler package, Flower, or sff2fastq.
* You can also manually map the reads yourself.  This could be a good option, but it could introduce some bias from the mapper software you use.  If you make the bam file yourself, place it into the `project/bam` directory and name it exactly the way that Lyve-SET would expect it.  In this way, you trick Lyve-SET into thinking that it already produced the bam file so that it will not recreate it and will include it in downstream analysis.
* You can also assemble the genome and include the assembly only.  See below on how to include an assembly.  This could potentially introduce an assembler bias.

##How do I include an assembly in the analysis?

Copy or symlink the assembly into the `project/asm` folder.  Lyve-SET uses the `samtools` program `wgsim` to simulate reads without introducing any errors.  These simulated reads will appear as `wgsim.fastq.gz` files in the reads directory.  Unfortunately, including assembly files in this way could introduce assembly-based errors.

High-quality-ness
=================

What are the different ways that SNPs in Lyve-SET have high confidence?
-----------------------------------------------------------------------
* **Detection of troublesome regions** such that they are not considered in hqSNP analysis.  Currently in v1.0, only phage genes are detected; however other databases could be added in the future, and also I am open to other suggestions.  Users can also specify a BED-formatted file to describe regions to mask.
* Only **unambiguous mapping** allowed
* Default **75% consensus** and **4x** coverage thresholds.  These options can be changed when you launch Lyve-SET.  Certain presets (found in `config/presets.conf`) can alter these thresholds.  For example, `--preset salmonella_enterica` invokes 20x/95%.  Additionally, each SNP must have at least two supporting reads in each direction.
* Mechanisms to **remove clustered SNPs** -- not on by default however.
* **Maximum likelihood** phylogeny reconstruction. Ascertainment bias is also considered through RAxML v8.
