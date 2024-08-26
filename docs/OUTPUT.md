Below is a visualization of the workflow with output files.
Then in the next section, a table with the description of all output files.

Visualization of output files
=============================

```mermaid
flowchart TD
  subgraph LOGDIR ["log directory"]
    LOG["main log file"]
    subgraph LOGSUBDIR ["SGELK subfolder"]
      LOGS["launch_set.pl.$$.log: individual log files"]
      JOBS["qsub.$$.pl: individual jobs to execute"]
      FINISHED["launch_set.pl.$$.finished: individual jobs that were finished"]
      SUBMITTED["launch_set.pl.$$.submitted: individual jobs that were submitted"]
      RUNNING["launch_set.pl.$$.running: individual jobs that are running"]
    end
  end
  
  SET_MANAGE_CREATE[--create] --> |Create all directories such as reads, reference| READSDIR
  SET_MANAGE_READS[--add-reads] --> |Create all directories such as reads, reference| READSDIR
  SET_MANAGE_ASM[--add-assembly] --> |Add an assembly| ASMDIR
  SET_MANAGE_REF[--change-reference] --> |Add a reference genome| REFDIR
  SET_MANAGE --- SET_MANAGE_CREATE
  SET_MANAGE --- SET_MANAGE_READS
  SET_MANAGE --- SET_MANAGE_ASM
  SET_MANAGE --- SET_MANAGE_REF
  subgraph READSDIR ["reads directory"]
    direction LR;
    SEQUENCER["R1 R2 fastq files"]
    RAW["Raw reads"]
    CLEANED["Cleaned reads"]

    RAW --> |run_assembly_trimClean.pl| CLEANED
    SEQUENCER --> |shuffleSplitReads.pl| RAW

  end
  subgraph REFDIR ["reference directory"]
    direction LR;
    REF["Reference genome"]
    UNMASKEDBED["unmaskedRegions.bed"]
    MASKEDBED["maskedRegions.bed"]

    REF --> |findPhages.pl| MASKEDBED
    MASKEDBED --> |invert| UNMASKEDBED
  end
  subgraph ASMDIR ["asm directory"]
    INASM["input assemblies"]
  end
  LAUNCH_SMALT{{launch_smalt.pl}}

  ASMDIR --> |samtools wgsim| READSDIR
  READSDIR --> |launch_smalt.pl| LAUNCH_SMALT
  REFDIR --> |launch_smalt.pl; only accept unmasked regions| LAUNCH_SMALT
  LAUNCH_SMALT --> |make a bam, once per genome| BAMDIR
  subgraph BAMDIR ["bam directory"]
    direction LR;
    BAM["*.bam"]
    BAMIDX["*.bam.bai"]
    BAMCLIFFS["*.bam.cliffs.bed"]

    BAM --> |set_findCliffs.pl| BAMCLIFFS
    BAM --> |samtools index| BAMIDX
  end
  LAUNCH_VARSCAN{{launch_varscan.pl}}

  BAMDIR --> |launch_varscan.pl \nexclude any regions in *.bam.cliffs.bed\nAccept SNPs at %consenus, X depth, fwd/rev support| LAUNCH_VARSCAN
  REFDIR --> |launch_varscan.pl\ndo not reprocess with unmaskedRegions.bed \nb/c launch_smalt.pl already used it| LAUNCH_VARSCAN
  LAUNCH_VARSCAN --> |make a vcf, once per genome| VCFDIR
  subgraph VCFDIR ["VCF directory"]
    direction LR;
    VCF["vcf files"]
    VCFIDX["vcf index files"]

    VCF --> |bcftools index| VCFIDX
  end
  VCFDIR --> |set_mergeVcf.sh: combine all Vcfs into \nout.pooled.vcf.gz and out.pooled.snps.vcf.gz| MSADIR
  subgraph MSADIR ["MSA directory"]
    direction LR;
    MERGEVCF{{set_mergeVcf.sh}}
    VCFPOOLED["out.pooled.vcf.gz"]
    VCFSNPSPOOLED["out.pooled.snps.vcf.gz"]
    SNPMATRIX["out.snpmatrix.tsv"]
    FILTEREDMATRIX["out.filteredMatrix.tsv"]
    FULLALN["out.aln.fasta"]
    INFORMATIVEALN["out.informative.fasta"]
    TREE["out.RAxML_bipartitions; tree.dnd"]
    PAIRWISE["out.pairwise.tsv"]
    PAIRWISEMATRIX["out.pairwiseMatrix.tsv"]

    MERGEVCF --> VCFPOOLED
    MERGEVCF --> VCFSNPSPOOLED
    VCFPOOLED --> |pooledToMatrix.sh| SNPMATRIX
    SNPMATRIX --> |filterMatrix.pl| FILTEREDMATRIX
    SNPMATRIX --> |matrixToAlignment.pl| FULLALN
    FILTEREDMATRIX --> |matrixToAlignment.pl| INFORMATIVEALN
    FULLALN --> |pairwiseDistances.pl| PAIRWISE
    PAIRWISE --> |pairwiseTo2d.pl| PAIRWISEMATRIX
    INFORMATIVEALN --> |launch_raxml.sh| TREE
  end
```

Output files
============
| File    |   Description  | Notes |
|:--------|:---------------|:------|
|project/msa | The multiple sequence alignment directory | Most of the output files you want are here like the multiple sequence alignment and the phylogeny|
|`project/msa/out.pooled.vcf.gz` | The pooled VCF file created from `bcftools merge` | 
|`project/msa/out.pooled.snps.vcf.gz` | SNPs vcf | The same data as `out.pooled.vcf.gz` but filtered to SNPs only. |
|`project/msa/out.pooled.vcf.gz.tbi`, `out.pooled.snps.vcf.gz.tbi` | the tabix index file for each VCF | 
|`project/msa/out.snpmatrix.tsv` | The `bcftools query` output | This file is essentially the main SNP matrix and describes the position and allele for each genome.  Each allele is in the genotype (GT) format, as specified in the vcf format specification |
|`project/msa/out.filteredMatrix.tsv` | The filtered `bcftools query` output | After `out.snpmatrix.tsv` is generated, this file describes remaining SNPs after some are filtered out, usually because the `--allowedFlanking` option in `launch_set.pl`, `--allowed` in `filterMatrix.pl`, or similar parameters in other scripts |
|`project/msa/out.aln.fasta` | The output alignment file in fasta format. | Make any changes to this file before running a phylogeny |program.  Do not use `out.informative.fasta` to make edits because positions might come and go and therefore you might lose resolution. After any edits, use `removeUninformativeSites.pl` to re-create `out.informative.fasta`  |
| `project/msa/out.informative.fasta` | The alignment after removing uninformative columns (ambiguities, invariants, gaps) | Do not make any changes to this file before running a phylogeny. Make the changes in `out.aln.fasta` |
| `project/msa/out.RAxML_bipartitions` | RAxML-generated tree in newick format | 
| `project/msa/tree.dnd` | Symlink to `out.RAxML_bipartitions`| 
| `project/msa/out.pairwise.tsv` | Pairwise distances file | Format: tab-delimited with three columns: genome1, genome2, hqSNP distance |
| `project/msa/out.pairwiseMatrix.tsv` | Pairwise distances matrix | The same data as `out.pairwise.tsv`, but in a 2-d matrix. Generated with `pairwiseTo2d.pl`. |
|project/log| Log files|
|`project/log/launch_set.log`    | The main log file |
|project/asm, project/reads || The input assemblies and reads. |
|project/reference | Where the reference fasta file is|
|`project/reference/maskedRegions.bed` | Regions of the reference genome that is masked for analysis. |
|project/reference/maskedRegions | BED-formatted files that describe regions that should be masked in the reference genome.|  You may also create your own file that can have any filename with extension `.bed`. This file can describe your manually-chosen regions that should be masked.  These regions will be incorporated into `project/reference/maskedRegions.bed`.|
|`project/reference/maskedRegions/phages.bed`| BED-formatted file describing predicted phage sites||
|project/bam| Output bam files are here|
|`project/bam/*.sorted.bam` | Sorted bam files | The query and reference name are encoded in the filename; many times the reference name will just be called "reference." |
|`project/bam/*.sorted.bam.bai` | Samtools index file |
|`project/bam/*.sorted.bam.cliffs.bed` | Files describing genome depth cliffs | These are only present if you specified `--mask-cliffs` |
|project/vcf |VCF files|Have the same file format as the `*.sorted.bam` files, so that they can be matched easily when running Lyve-SET. These files are sorted with vcftools and compressed with bgzip.|
|`project/vcf/*.vcf.gz`|VCF files ||
|`project/vcf/*.vcf.gz.tbi`| Tabix index files|
