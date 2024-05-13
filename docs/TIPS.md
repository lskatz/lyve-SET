Tips and Tricks
===============

Here are just some tips and tricks that I have used or that others have contributed

Masking a region in your reference genome
-----------------------------------------
Yes, there is actually a mechanism to manually mask troublesome regions in the reference genome!  Under `project/reference/maskedRegions`, create a file with an extension `.bed`.  This file has at least three columns: `contig`, `start`, `stop`.  BED is a standard file format and is better described here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1 

The Lyve-SET phage-finding tool that uses PHAST actually puts a `phages.bed` file into that directory.  In the course of the pipeline, Lyve-SET will use any BED files in that directory to 1) ignore any reads that are mapped entirely in those regions and 2) ignore any SNPs that are found in those regions.  In the future, Lyve-SET will also use individualized BED files in the bam directory to mask SNPs found on a per-genome basis.

Using multiple processors on a single-node machine
--------------------------------------------------

Unfortunately if you are not on a cluster, then Lyve-SET will only work on a single node.  Here are some ways of using xargs to speed up certain steps of Lyve-SET. Incorporating these changes or something similar is on my todo list but for now it is easier to post them.

### Making pileups

In the case where you want to generate all the pileups on one node using SNAP using 24 cores

    ls reads/*.fastq.gz | xargs -n 1 -I {} basename {} | xargs -P 24 -n 1 -I {} launch_smalt.pl -f reads/{} -b bam/{}-2010EL-1786.sorted.bam -t tmp -r reference/2010EL-1786.fasta --numcpus 1

### Calling SNPs

In the case where you have all pileups finished and want to call SNPs on a single node.  This example uses 24 cpus.  At the end of this example, you will still need to sort the VCF files (`vcf-sort`), compress them with `bgzip`, and index them with `tabix`.

    # Call SNPs into vcf file
    ls bam/*.sorted.bam | xargs -n 1 -I {} basename {} .sorted.bam | xargs -P 24 -n 1 -I {} sh -c "/home/lkatz/bin/Lyve-SET/scripts/launch_varscan.pl bam/{}.sorted.bam --tempdir tmp --reference reference/2010EL-1786.fasta --altfreq 0.75 --coverage 10 > vcf/{}.vcf"
    # sort/compress/index
    cd vcf; 
    ls *.vcf| xargs -I {} -P 24 -n 1 sh -c "vcf-sort < {} > {}.tmp && mv -v {}.tmp {} && bgzip {} && tabix {}.gz"

## SNP interrogation

### What are the differences between just two isolates?

Sometimes you just want to know the SNPs that define one isolate, or maybe find which SNPs define which clades.
Here is a method to show unique SNPs between two isolates.
You can expand this idea to more isolates or otherwise refine it as needed.

```bash
$ cut -f 1,2,5,6 out.snpmatrix.tsv | awk '$3!=$4 && $3!="N" && $4!="N"' | head -n 3
# [1]CHROM      [2]POS  [5]sample1:GT   [6]sample2:GT
gi|9626243|ref|NC_001416.1|     403     A       G
gi|9626243|ref|NC_001416.1|     753     A       G
```

In this example, only the first two sites are shown (`head -n 3`)
and show that at positions 403 and 753, there are differences between sample1 and sample2.
Sites with `N` are ignored with logic like `$3!="N"`.
Otherwise `awk` is just selecting for rows where the samples differ (`$3!=$4`)
`cut` is used to grab only the position columns 1 and 2 and then to keep only two samples at columns 5 and 6.

### What are the differences between two groups of isolates?

To do this, we can run `bcftools isec` in combination with `bcftools filter` and `bcftools view`.
In this section, we can use the lambda subset to pretend that there are two clades with clade1 having sample1 and sample2,
and clade2 having samples 3 and 4.

Go into the `msa` subfolder with `cd msa`.

    # Make files with sample1 and sample2, and then sample3 and sample4
    echo -e "sample1\nsample2" > clade1.txt
    echo -e "sample3\nsample4" > clade2.txt

Use the file that already contains SNPs for all samples and divide it into the two clades

    bcftools filter -i 'ALT!="N"' out.pooled.snps.vcf.gz | \
      bcftools view -S clade1.txt | \
      bgzip -c > clade1.snps.vcf.gz
    bcftools filter -i 'ALT!="N"' out.pooled.snps.vcf.gz | \
      bcftools view -S clade2.txt | \
      bgzip -c > clade2.snps.vcf.gz

Index the new VCF files

    tabix clade1.snps.vcf.gz
    tabix clade2.snps.vcf.gz

Find the intersection of all SNPs between the two clades

    bcftools isec -p isec clade1.snps.vcf.gz clade2.snps.vcf.gz

Go into the new `isec` folder and format all the vcf files

    cd isec
    for i in *.vcf; do
      bgzip $i
      tabix $i.gz
    done

View the notes on which file is which. Typically, `0000.vcf.gz` is going to be sites exclusive to only clade1.
`0001.vcf.gz` is going to be sites exclusive to only clade2.
`0002.vcf.gz` is records in clade1 that are also in clade2.
`0003.vcf.gz` is records in clade2 that are also in clade1.

    cat README.txt

    This file was produced by vcfisec.
    The command line was:   bcftools isec  -p isec clade1.snps.vcf.gz clade2.snps.vcf.gz
    
    Using the following file names:
    isec/0000.vcf   for records private to  clade1.snps.vcf.gz
    isec/0001.vcf   for records private to  clade2.snps.vcf.gz
    isec/0002.vcf   for records from clade1.snps.vcf.gz shared by both      clade1.snps.vcf.gz clade2.snps.vcf.gz
    isec/0003.vcf   for records from clade2.snps.vcf.gz shared by both      clade1.snps.vcf.gz clade2.snps.vcf.gz

Merge back into one happy VCF file

    bcftools merge 0002.vcf.gz 0003.vcf.gz > merged.vcf
    bgzip merged.vcf
    tabix merged.vcf.gz

View it. These are now SNPs that _could_ differentiate clade1 from clade2, but not necessarily 100%.
For example, a SNP might only occur in sample1 but not sample2 even though it's the same clade.

    query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' merged.vcf.gz | column -ts $'\t' | less -S

In our example with lambda, with its star phylogeny, however, we do not have SNPs that are 100%.
Hopefully your results differ in this regard :) 

    # [1]CHROM                   [2]POS  [3]REF  [4]ALT  [5]sample1:GT  [6]sample2:GT  [7]sample3:GT  [8]sample4:GT
    gi|9626243|ref|NC_001416.1|  403     G       A       1/1            0/0            0/0            0/0
    gi|9626243|ref|NC_001416.1|  550     G       A       0/0            0/0            0/0            1/1
    gi|9626243|ref|NC_001416.1|  586     C       G       0/0            0/0            0/0            1/1
    gi|9626243|ref|NC_001416.1|  753     G       A       1/1            0/0            0/0            0/0
    gi|9626243|ref|NC_001416.1|  1019    C       T       0/0            0/0            0/0            1/1

### SNP counting

How many sites are in the Lyve-SET analysis?  Count the number of lines in `out.snpmatrix.tsv`. Subtract 1 for the header.

How many SNPs are in the Lyve-SET analysis? Count the number of lines in `out.filteredMatrix.tsv`. Subtract 1 for the
header.

I want to count the number of sites as defined as X.  Use `out.snpmatrix.tsv` and `filterMatrix.pl` to filter it your way. If you are an advanced user, you can use `bcftools query` on `out.pooled.vcf.gz`.

### Why did a SNP fail?

To see why any SNP failed, view the vcf files in the VCF directory.  Parse them with `bcftools`.

#### View the FILTER column

    $ bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n' sample1.fastq.gz-reference.vcf.gz | head
    gi|9626243|ref|NC_001416.1|     1       G       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     2       G       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     3       G       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     4       C       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     5       G       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     6       G       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     7       C       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     8       G       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     9       A       N       DP10;AF0.75
    gi|9626243|ref|NC_001416.1|     10      C       N       DP10;AF0.75

#### View MORE information with bcftools

    $ bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%TGT\t%DP]\n' sample1.fastq.gz-reference.vcf.gz | head
    gi|9626243|ref|NC_001416.1|     1       G       N       DP10;AF0.75     N/N     2
    gi|9626243|ref|NC_001416.1|     2       G       N       DP10;AF0.75     N/N     4
    gi|9626243|ref|NC_001416.1|     3       G       N       DP10;AF0.75     N/N     4
    gi|9626243|ref|NC_001416.1|     4       C       N       DP10;AF0.75     N/N     4
    gi|9626243|ref|NC_001416.1|     5       G       N       DP10;AF0.75     N/N     4
    gi|9626243|ref|NC_001416.1|     6       G       N       DP10;AF0.75     N/N     5
    gi|9626243|ref|NC_001416.1|     7       C       N       DP10;AF0.75     N/N     5
    gi|9626243|ref|NC_001416.1|     8       G       N       DP10;AF0.75     N/N     6
    gi|9626243|ref|NC_001416.1|     9       A       N       DP10;AF0.75     N/N     6
    gi|9626243|ref|NC_001416.1|     10      C       N       DP10;AF0.75     N/N     6

#### View the definitions for FILTER codes

For example, `DP10` indicates that the user set 10x as a threshold and the site has less than 10x coverage.
    
    $ zgrep "##FILTER" sample1.fastq.gz-reference.vcf.gz
    ##FILTER=<ID=DP10,Description="Depth is less than 10, the user-set coverage threshold">
    ##FILTER=<ID=RF0.75,Description="Reference variant consensus is less than 0.75, the user-set threshold">
    ##FILTER=<ID=AF0.75,Description="Allele variant consensus is less than 0.75, the user-set threshold">
    ##FILTER=<ID=isIndel,Description="Indels are not used for analysis in Lyve-SET">
    ##FILTER=<ID=masked,Description="This site was masked using a bed file or other means">
    ##FILTER=<ID=str10,Description="Less than 10% or more than 90% of variant supporting reads on one strand">
    ##FILTER=<ID=indelError,Description="Likely artifact due to indel reads at this position">

### Find constant sites

Lyve-SET focuses on variable sites, but what if you wanted to look at how many constant sites there are?
This is at least useful for when you need to tell BEAST or other similar programs what your constant sites are.

To get the full alignment, you can run the following `bcftools` command, 
followed by some parsing

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t[%TGT\t]\n' --print-header out.pooled.vcf.gz > full.tsv.unrefined.tmp

# Get the first nucleotide of every nucleotide call to convert from diploid to haploid
perl -lane '
  BEGIN{
  # print the header
  $line=<>;
  chomp($line);
  print $line;
  $numFields=scalar(split(/\t/,$line));
  $lastIndex=$numFields-1;
}
for(@F[3..$lastIndex]){
  $_=substr($_,0,1);
}
print join("\t",@F);
' < full.tsv.unrefined.tmp > full.tsv
# Now you have full.tsv

# Parse full.tsv to get counts of constant sites
tail -n +2 full.tsv | perl -MData::Dumper -lane '
  # Print a dot every 100k lines
  print STDERR "Looked at $. lines" if($. % 100000 == 0);

  # Remove the contig, pos, ref fields
  splice(@F,0,3);
  # Pretend the first sample is the reference
  # and remove it from the list of samples with shift()
  $ref=shift(@F);
  # Skip if "reference" base is ambiguous
  next if($ref eq "N");
  # The site is constant until proven otherwise
  $is_constant=1;
  # Loop through all samples after the first
  for my $nt(@F){
    # If the samples nucleotide is not equal to the first samples,
    # then label as not constant and go to the next site
    if($nt ne $ref){
      $is_constant=0;
      last;
    }
  }
  $const{$ref} += $is_constant;
  END{
    print Dumper \%const;
    print STDERR "\n";
  }
'
```

This will give you output similar to

```text
$VAR1 = {
          'A' => 638550,
          'T' => 683152,
          'C' => 443617,
          '.' => 0,
          'G' => 412973
        };
```

Where you get counts for every nucleotide where it was constant at any given site in the set of samples.
You might be able to put it into BEAST with some syntax similar to (but not necessarily):

```xml
  <alignment dataType="nucleotide">
    <sequence idref="taxon1" value="G" count="412973"/>
    <sequence idref="taxon1" value="C" count="443617"/>
    <sequence idref="taxon1" value="T" count="683152"/>
    <sequence idref="taxon1" value="." count="0"/>
    <sequence idref="taxon1" value="A" count="638550"/>
    <!-- Add more sequence elements for each taxon and its corresponding nucleotide with their respective counts -->
  </alignment>
```

And then, lastly, to get the full alignment including constant sites, you can transform `full.tsv` like so

```bash
matrixToAlignment.pl full.tsv > full.fasta
```

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
    
