Tips and Tricks
===============

Here are just some tips and tricks that I've used or that others have contributed

Using multiple processors on a single-node machine
--------------------------------------------------

Unfortunately if you are not on a cluster, then Lyve-SET will only work on a single node.  Here are some ways of using xargs to speed up certain steps of Lyve-SET. Incorporating these changes or something similar is on my todo list but for now it is easier to post them.

###Making pileups

In the case where you want to generate all the pileups on one node using SNAP using 24 cores

    ls reads/*.fastq.gz | xargs -n 1 -I {} basename {} | xargs -P 24 -n 1 -I {} launch_snap.pl -f reads/{} -b bam/{}-2010EL-1786.sorted.bam -t tmp -r reference/2010EL-1786.fasta --numcpus 1

###Calling SNPs

In the case where you have all pileups finished and want to call SNPs on a single node.  This example uses 24 cpus.  At the end of this example, you'll still need to sort the VCF files (`vcf-sort`), compress them with `bgzip`, and index them with `tabix`.

    ls bam/*.sorted.bam | xargs -n 1 -I {} basename {} .sorted.bam | xargs -P 24 -n 1 -I {} sh -c "/home/lkatz/bin/Lyve-SET/scripts/launch_varscan.pl bam/{}.sorted.bam --tempdir tmp --reference reference/2010EL-1786.fasta --altfreq 0.75 --coverage 10 > vcf/{}.vcf"
