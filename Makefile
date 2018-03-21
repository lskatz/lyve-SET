# Author: Lee Katz <lkatz@cdc.gov>
# Lyve-SET
 
PREFIX  := $(PWD)
PROFILE := $(HOME)/.bashrc
PROJECT := "setTestProject"
NUMCPUS := 1
SHELL   := /bin/bash

# Derived variables
TMPDIR := build
TMPTARFILE=$(TMPDIR)/$(TARFILE)

# Style variables
T= "	"
T2=$(T)$(T)

###################################

default: install 

help:
	@echo "Please see README.md for additional help"

all: install env

install: install-prerequisites
	@echo "Don't forget to set up update PATH to $(PREFIX)/scripts";
	@echo "'make env' performs this step for you"
	@echo "DONE: installation of Lyve-SET complete."

install-prerequisites: install-mkdir install-vcftools install-CGP install-SGELK install-varscan install-phast install-samtools install-bcftools install-smalt install-snap install-bwa install-raxml install-perlModules install-config install-snpEff install-eutils
	@echo DONE installing prerequisites

clean: clean-tmp clean-symlinks clean-vcftools clean-CGP clean-SGELK clean-varscan clean-phast clean-samtools clean-bcftools clean-smalt clean-bwa clean-snap clean-raxml clean-perlModules clean-config clean-snpEff clean-eutils
	@echo "Remember to remove the line with PATH and Lyve-SET from $(PROFILE)"

install-mkdir:
	-mkdir build lib scripts

clean-tmp:
	rm -rfv $(TMPDIR)/*

clean-symlinks:
	find scripts -maxdepth 1 -type l -exec rm -vf {} \;

install-SGELK:
	git clone https://github.com/lskatz/Schedule--SGELK.git $(TMPDIR)/Schedule
	-mkdir -p lib/Schedule
	mv -v $(TMPDIR)/Schedule/SGELK.pm lib/Schedule/
	mv -v $(TMPDIR)/Schedule/README.md lib/Schedule/
	mv -v $(TMPDIR)/Schedule/.git lib/Schedule/

clean-SGELK:
	rm -rfv lib/Schedule

install-CGP:
	# CGP scripts that are needed and that don't depend on CGP libraries
	git clone https://github.com/lskatz/cg-pipeline lib/cg-pipeline
	ln -rs $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_isFastqPE.pl scripts/
	ln -rs $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_trimClean.pl scripts/
	ln -rs $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_shuffleReads.pl scripts/
	ln -rs $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_removeDuplicateReads.pl scripts/
	ln -rs $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_readMetrics.pl scripts/
	ln -rs $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_metrics.pl scripts/

clean-CGP:
	rm -rvf lib/cg-pipeline
	rm -vf scripts/run_assembly_*

install-vcftools:
	wget --max-redirect 50 'http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F&ts=1409260024&use_mirror=ufpr' -O $(TMPDIR)/vcftools_0.1.12b.tar.gz
	cd $(TMPDIR) && \
	tar zxvf vcftools_0.1.12b.tar.gz
	mv $(TMPDIR)/vcftools_0.1.12b lib/
	cd lib/vcftools_0.1.12b &&\
    $(MAKE) --directory=. MAKEFLAGS=""
	ln -s $(PREFIX)/lib/vcftools_0.1.12b/perl/vcf-sort scripts/
	ln -s $(PREFIX)/lib/vcftools_0.1.12b/perl/Vcf.pm lib/

clean-vcftools:
	rm -rvf lib/vcftools_0.1.12b
	rm -vf scripts/vcf-sort
	rm -vf lib/Vcf.pm


install-varscan:
	wget --max-redirect 50 'http://downloads.sourceforge.net/project/varscan/VarScan.v2.3.7.jar' -O lib/varscan.v2.3.7.jar

clean-varscan:
	rm -vf lib/varscan.v2.3.7.jar

install-snpEff:
	wget --max-redirect 50 http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip -O lib/snpEff_latest_core.zip
	cd lib && unzip -o snpEff_latest_core.zip
	mv lib/snpEff/snpEff.jar lib/
	mv lib/snpEff/snpEff.config config/snpEff.conf
	rm -rf lib/snpEff_latest_core.zip
	rm -rf lib/snpEff

clean-snpEff:
	rm -rvf lib/snpEff.jar

install-phast: check-blast
	mkdir -p lib/phast
	wget http://phast.wishartlab.com/phage_finder/DB/prophage_virus.db -O lib/phast/phast.faa
	makeblastdb -in lib/phast/phast.faa -dbtype prot

clean-phast:
	rm -rvf lib/phast

install-samtools:
	wget 'https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2' -O $(TMPDIR)/samtools-1.3.1.tar.bz2
	wget 'https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2' -O $(TMPDIR)/htslib-1.3.2.tar.bz2
	cd $(TMPDIR) && tar jxvf htslib-1.3.2.tar.bz2 && tar jxvf samtools-1.3.1.tar.bz2
	mv $(TMPDIR)/htslib-1.3.2 $(TMPDIR)/samtools-1.3.1 lib
	cd lib/htslib-1.3.2 && $(MAKE)
	cd lib/samtools-1.3.1 && $(MAKE) DFLAGS="-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_CURSES_LIB=0" LIBCURSES="" 
	ln -rsf $(PREFIX)/lib/samtools-1.3.1/samtools scripts/
	ln -rsf $(PREFIX)/lib/samtools-1.3.1/misc/wgsim scripts/
	ln -rsf $(PREFIX)/lib/htslib-1.3.2/bgzip scripts
	ln -rsf $(PREFIX)/lib/htslib-1.3.2/tabix scripts

clean-samtools:
	rm -rvf lib/samtools*
	rm -rvf lib/htslib*
	rm -vf scripts/{samtools,wgsim,bgzip,tabix}

install-bcftools:
	wget 'https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2' -O $(TMPDIR)/bcftools-1.3.1.tar.bz2
	cd $(TMPDIR) && tar jxvf bcftools-1.3.1.tar.bz2
	mv $(TMPDIR)/bcftools-1.3.1 lib/bcftools-1.3.1
	cd lib/bcftools-1.3.1 && $(MAKE)
	ln -rs $(PREFIX)/lib/bcftools-1.3.1/bcftools scripts
	ln -rs $(PREFIX)/lib/bcftools-1.3.1/vcfutils.pl scripts

clean-bcftools:
	rm -rfv $(PREFIX)/lib/bcftools*
	rm -vf $(PREFIX)/scripts/vcfutils.pl scripts/bcftools

install-smalt:
	wget --max-redirect 50 --continue 'https://downloads.sourceforge.net/project/smalt/smalt-0.7.6-static.tar.gz' -O $(TMPDIR)/smalt-0.7.6-static.tar.gz
	cd $(TMPDIR) && tar zxvf smalt-0.7.6-static.tar.gz
	mv $(TMPDIR)/smalt-0.7.6 lib/
	cd lib/smalt-0.7.6 && ./configure --prefix $(PREFIX)/lib/smalt-0.7.6
	$(MAKE) --directory lib/smalt-0.7.6
	$(MAKE) --directory lib/smalt-0.7.6 install
	ln -rsv $(PREFIX)/lib/smalt-0.7.6/bin/smalt scripts/

clean-smalt:
	rm -rvf lib/smalt*
	rm -vf scripts/smalt

install-bwa:
	wget --max-redirect 50 --continue 'https://astuteinternet.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2' -O $(TMPDIR)/bwa-0.7.17.tar.bz2
	cd $(TMPDIR) && tar jxvf bwa-0.7.17.tar.bz2
	mv $(TMPDIR)/bwa-0.7.17 lib/
	$(MAKE) --directory lib/bwa-0.7.17
	ln -rsv $(PREFIX)/lib/bwa-0.7.17/bwa scripts/

clean-bwa:
	rm -rvf lib/bwa*
	rm -vf scripts/bwa

install-snap:
	git clone https://github.com/amplab/snap.git lib/snap -b v1.0beta.18
	cd lib/snap &&  $(MAKE)
	cp -vn lib/snap/snap-aligner  scripts/snap
	cd lib/snap && $(MAKE) clean && $(MAKE) CXXFLAGS="-DLONG_READS -O3 -Wno-format -MMD -ISNAPLib -msse"
	cp -vn lib/snap/snap-aligner  scripts/snapxl

clean-snap:
	rm -rvf lib/snap scripts/snap scripts/snapxl

install-raxml: lib/standard-RAxML-8.1.16/raxmlHPC lib/standard-RAxML-8.1.16/raxmlHPC-PTHREADS

lib/standard-RAxML-8.1.16/raxmlHPC:
	wget 'https://github.com/stamatak/standard-RAxML/archive/v8.1.16.tar.gz' -O $(TMPDIR)/raxml_v8.1.16.tar.gz
	cd $(TMPDIR) && tar zxvf raxml_v8.1.16.tar.gz
	mv $(TMPDIR)/standard-RAxML-8.1.16 lib
	$(MAKE) --directory lib/standard-RAxML-8.1.16 -f Makefile.gcc
	ln -rsv $(PREFIX)/$@ scripts/
lib/standard-RAxML-8.1.16/raxmlHPC-PTHREADS: lib/standard-RAxML-8.1.16/raxmlHPC
	$(MAKE) --directory lib/standard-RAxML-8.1.16 -f Makefile.PTHREADS.gcc
	ln -rsv $(PREFIX)/$@ scripts/

clean-raxml:
	rm -rvf lib/standard-RAxML-8.1.16
	rm -vf  scripts/raxmlHPC scripts/raxmlHPC-PTHREADS

install-edirect:
	cd build && perl -MNet::FTP -e '$$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1); $$ftp->login; $$ftp->binary; $$ftp->get("/entrez/entrezdirect/edirect.zip");' && unzip -u -q edirect.zip && rm edirect.zip
	cd build/edirect && sh setup.sh
	mv -v build/edirect lib

clean-edirect:
	rm -rvf lib/edirect

install-perlModules:
	@echo "Installing Perl modules using cpanminus"
	#for package in Config::Simple File::Slurp Math::Round Number::Range Statistics::Distributions Statistics::Descriptive Statistics::Basic Graph::Centrality::Pagerank String::Escape Statistics::LineFit; do
	for package in Config::Simple File::Slurp Math::Round Number::Range Statistics::Distributions Statistics::Basic Graph::Centrality::Pagerank String::Escape Statistics::LineFit Array::IntSpan; do \
	  perl scripts/cpanm --self-contained -L lib $$package; \
		if [ $$? -gt 0 ]; then exit 1; fi; \
	done;
	@echo "Done with Perl modules"

clean-perlModules:
	@echo "Perl modules were installed using CPAN which doesn't have an uninstalling mechanism"

install-config:
	cp -vn config/original/*.conf config/

clean-config:
	rm -vf config/*.conf

env:
	echo "#Lyve-SET" >> $(PROFILE)
	echo "export PATH=\$$PATH:$(PREFIX)/scripts" >> $(PROFILE)
	echo "export PERL5LIB=\$$PERL5LIB:$(PREFIX)/lib" >> $(PROFILE)

# https://www.ncbi.nlm.nih.gov/books/NBK179288/
install-eutils:
	wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.zip -O build/edirect.zip
	unzip -d build build/edirect.zip
	mv build/edirect lib/edirect
	rm build/edirect.zip
	lib/edirect/setup.sh

clean-eutils:
	rm -rvf lib/edirect

test:
	set_test.pl lambda lambda --numcpus $(NUMCPUS) 

check: check-sys check-Lyve-SET-PATH check-CGP-assembly check-Lyve-SET check-PERL check-smalt check-raxml check-phyml check-blast
	@echo --OK
check-sys:
	@F=$$(which cat) && echo "Found $$F"
	@F=$$(which gzip) && echo "Found $$F"
	@F=$$(which perl) && echo "Found $$F"
check-smalt:
	@F=$$(which smalt 2>/dev/null) && echo "Found smalt at $$F"
check-raxml:
	@ RAXML=$$((which raxml || which raxmlHPC-PTHREADS) 2>/dev/null) && echo "Found raxml at $$RAXML"
check-phyml:
	@ PHYML=$$((which phyml || which phyml_linux_64 ) 2>/dev/null) && echo "Found phyml at $$PHYML"
check-CGP-assembly:
	@which run_assembly_shuffleReads.pl run_assembly_trimClean.pl run_assembly_isFastqPE.pl >/dev/null
check-Lyve-SET-PATH:
	@echo Checking that your path includes Lyve-SET/scripts
	@which launch_set.pl >/dev/null
check-Lyve-SET:
	@echo Checking that the SET executables are present
	@export PATH=$$PATH:$(PREFIX)/scripts && which launch_set.pl >/dev/null
check-PERL:
	@echo Checking for perl multithreading
	@perl -Mthreads -e 1
	@echo Checking for perl modules
	@echo "Looking for File::Slurp"
	@perl -I $(PREFIX)/lib -MFile::Slurp -e 1
	@echo "Looking for String::Escape"
	@perl -I $(PREFIX)/lib -MString::Escape -e 1
check-blast:
	which blastx
	which blastp
	which makeblastdb
