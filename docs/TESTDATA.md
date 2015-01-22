To make a test dataset, create a directory with all lowercase with the the files SRA and NUCLEOTIDE.  These files are described below; you do not have to create any other files in the test data directory.

All directories must start off by containing the files SRA and NUCLEOTIDE which have identifiers for raw reads and for genome assemblies (draft or finished). The first identifier in NUCLEOTIDE will be the reference genome.

The final files/directories in each must be the following
asm/  NUCLEOTIDE  reads/  reference/  SRA
* asm/ has all the genome assemblies
* NUCLEOTIDE has all assembly identifiers. The reference genome is the first listing
* reads/ has all raw read fastq files
* reference/ contains a single reference genome assembly
* SRA has a listing of all raw read identifiers that should be downloaded

CREDITS
* Test lambda data were obtained from the CFSAN SNP pipeline at https://github.com/CFSAN-Biostatistics/snp-pipeline
* Other datasets were created between EDLB and CFSAN
