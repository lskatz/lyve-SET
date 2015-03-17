Installation
============

Requirements
------------
* **Must-have** and not installed with `make install`
  * **Perl, multithreaded**
  * **BLAST+**
  * **GIT**, **SVN** (for installation and updating)
* **Requirements installed with `make install`**
  * CG-Pipeline
  * Schedule::SGELK
  * VCF Tools and Vcf.pm
  * Samtools v1.1 or newer
    * bcftools
    * wgsim if you are simulating reads
    * bgzip
    * tabix
  * PHAST
  * Varscan
  * SNAP
  * RAxML
  * Smalt
* Optional "requirements"
  * PhyML
  * FreeBayes

Installation
------------
* `make install`
* `make check` - check and see if you have all the prerequisites
* `make test` - run a test phage dataset provided by CFSAN
* `make help` - for other `make` options
* `make clean` - clean up an old installation in preparation for a new installation
* `make install-*` - Many other installation options are available including but not limited to:
  * `make install-smalt`
  * `make install-CGP`
  * `make install-samtools`
* `make clean-*` - Every `make install` command comes with a `make clean` command, e.g.:
  * `make clean-CGP`

