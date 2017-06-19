# Visualization

A large question is, after you are finished with a Lyve-SET run, how do you visualize the results?  One of the advantages of Lyve-SET is its use of standard file formats.  Therefore most files can be visualized in standard software.  All output files are documented under [OUTPUT.md](OUTPUT.md).

## Tree

You can visualize the tree in many different tree drawing programs out there.  Some of my personal favorites are [MEGA](http://www.megasoftware.net) and [Figtree](http://tree.bio.ed.ac.uk/software/figtree).  When prompted, open `tree.dnd`.  Some commercial software includes BioNumerics, CLC, and Geneious.

## SNP positions

SNP positions are encoded in the vcf files under the vcf directory.  You can view them in any standard viewer such as [IGV](http://software.broadinstitute.org/software/igv).  On the command line, although difficult to visualize, you can use `bcftools` which is redistributed in Lyve-SET.

To get started on IGV, first load the reference genome assembly. Second, load the bam and vcf files.  You will immediately be able to browse the genome with these bam and vcf tracks.  However for the advanced features, there is a learning curve (but it is worth it).

## Read alignments

Read alignments are encoded in bam files under the bam directory. These can be viewed with [IGV](http://software.broadinstitute.org/software/igv) and many commercial software packages such as CLC and Geneious.  Additionally if you are command-line-inclined, you can view them with `samtools tview` which is redistributed with Lyve-SET.

## SNP distances

SNP distances might be the easiest to visualize. These can be viewed in LibreOffice or Microsoft Excel.  Simply open `out.pairewiseMatrix.tsv` and turn on conditional formatting.  This shows the distance heatmap.
