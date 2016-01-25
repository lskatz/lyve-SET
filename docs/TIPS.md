Tips and Tricks
===============

Here are just some tips and tricks that I've used or that others have contributed

Masking a region in your reference genome
-----------------------------------------
Yes, there is actually a mechanism to manually mask troublesome regions in the reference genome!  Under `project/reference/maskedRegions`, create a file with an extension `.bed`.  This file has at least three columns: `contig`, `start`, `stop`.  BED is a standard file format and is better described here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1 

The Lyve-SET phage-finding tool that uses PHAST actually puts a `phages.bed` file into that directory.  In the course of the pipeline, Lyve-SET will use any BED files in that directory to 1) ignore any reads that are mapped entirely in those regions and 2) ignore any SNPs that are found in those regions.  In the future, Lyve-SET will also use individualized BED files in the bam directory to mask SNPs found on a per-genome basis.

Using multiple processors on a single-node machine
--------------------------------------------------

Unfortunately if you are not on a cluster, then Lyve-SET will only work on a single node.  Here are some ways of using xargs to speed up certain steps of Lyve-SET. Incorporating these changes or something similar is on my todo list but for now it is easier to post them.

###Making pileups

In the case where you want to generate all the pileups on one node using SNAP using 24 cores

    ls reads/*.fastq.gz | xargs -n 1 -I {} basename {} | xargs -P 24 -n 1 -I {} launch_snap.pl -f reads/{} -b bam/{}-2010EL-1786.sorted.bam -t tmp -r reference/2010EL-1786.fasta --numcpus 1

###Calling SNPs

In the case where you have all pileups finished and want to call SNPs on a single node.  This example uses 24 cpus.  At the end of this example, you'll still need to sort the VCF files (`vcf-sort`), compress them with `bgzip`, and index them with `tabix`.

    # Call SNPs into vcf file
    ls bam/*.sorted.bam | xargs -n 1 -I {} basename {} .sorted.bam | xargs -P 24 -n 1 -I {} sh -c "/home/lkatz/bin/Lyve-SET/scripts/launch_varscan.pl bam/{}.sorted.bam --tempdir tmp --reference reference/2010EL-1786.fasta --altfreq 0.75 --coverage 10 > vcf/{}.vcf"
    # sort/compress/index
    cd vcf; 
    ls *.vcf| xargs -I {} -P 24 -n 1 sh -c "vcf-sort < {} > {}.tmp && mv -v {}.tmp {} && bgzip {} && tabix {}.gz"

Other manual steps in Lyve-SET
-------------------------------------------------

### From a set of VCFs to finished results

    mergeVcf.sh -o msa/out.pooled.vcf.gz vcf/*.vcf.gz # get a pooled VCF
    cd msa
    pooledToMatrix.sh -o out.bcftoolsquery.tsv out.pooled.vcf.gz  # Create a readable matrix
    filterMatrix.pl --noambiguities --noinvariant  < out.bcftoolsquery.tsv > out.filteredbcftoolsquery.tsv # Filter out low-quality sites
    matrixToAlignment.pl < out.filteredbcftoolsquery.tsv > out.aln.fas  # Create an alignment in standard fasta format
    set_processMsa.pl --numcpus 12 --auto --force out.aln.fas # Run the next steps in this mini-pipeline

### Manual steps in `set_processMsa.pl`

Hopefully all these commands make sense but please tell me if I need to expound.

    cd msa
    removeUninformativeSites.pl --gaps-allowed --ambiguities-allowed out.aln.fas > /tmp/variantSites.fasta
    pairwiseDistances.pl --numcpus 12 /tmp/variantSites.fasta | sort -k3,3n | tee pairwise.tsv | pairwiseTo2d.pl > pairwise.matrix.tsv && rm /tmp/variantSites.fasta
    set_indexCase.pl pairwise.tsv | sort -k2,2nr > eigen.tsv # Figure out the most "connected" genome which is the most likely index case
    launch_raxml.sh -n 12 informative.aln.fas informative # create a tree with the suffix 'informative'
    applyFstToTree.pl --numcpus 12 -t RAxML_bipartitions.informative -p pairwise.tsv --outprefix fst --outputType averages > fst.avg.tsv  # look at the Fst for your tree (might result in an error for some trees, like polytomies)
    applyFstToTree.pl --numcpus 12 -t RAxML_bipartitions.informative -p pairwise.tsv --outprefix fst --outputType samples > fst.samples.tsv  # instead of average Fst values per tree node, shows you each repetition
    
## Interrogate the SNP calls

After everything is said and done, you might want to look for certain things like hotspots.  I'll start listing those tips and tricks here.  Because Lyve-SET depends on standard file formats, there are very standard ways to scrutinize the output files.

### Find pairwise SNPs using bcftools

    $ bcftools query --include '%TYPE="snp" && ALT!="N"' --print-header -s D7328,D7322 out.pooled.vcf.gz -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n'
