Reference genomes
=================

Can I have multiple reference genomes?
--------------------------------------
No.  However, you can add multiple assemblies to the asm directory.  Assemblies don't have the same error profile as reads and so you might expect some skew in the ultimate phylogeny.

How do I choose a reference genome?
-----------------------------------
The best reference genome for an outbreak is something in-clade.  You might even need to assemble the genome from your own reads before starting.  An outgroup is not as related to your clade by definition and so it is probably not the best genome to use.  The next best quality is a closed genome, or a genome with a high N50.

High-quality-ness
=================

What are the different ways that SNPs in Lyve-SET have high confidence?
-----------------------------------------------------------------------

High quality for SNPs indicates that the resulting phylogeny will be high-fidelity.  Although some SNPs are discarded that we are less sure about, the SNPs that we _are_ most sure about are retained, and the resulting phylogeny is the best inference.

* **Detection of troublesome regions** such that they are not considered in hqSNP analysis.  Currently in v1.0, only phage genes are detected; however other databases could be added in the future, and also I am open to other suggestions.  Users can also specify a BED-formatted file to describe regions to mask.
* Only **unambiguous mapping** allowed
* Default **75% consensus** and **10x** coverage thresholds.  These options can be changed when you launch Lyve-SET.
* Mechanisms to **remove clustered SNPs** -- not on by default however.
* **Maximum likelihood** phylogeny reconstruction. Ascertainment bias is also considered through RAxML v8.
