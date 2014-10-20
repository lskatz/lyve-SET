# Author: Lee Katz <lkatz@cdc.gov>
# Lyve-SET

PREFIX := $(PWD)
PROFILE := $(HOME)/.bashrc
VERSION := 0.9.2
PROJECT := "setTestProject"
NUMCPUS := 1
SHELL   := /bin/bash

# Derived variables
TMPDIR := $(PREFIX)/build
TMPTARFILE=$(TMPDIR)/$(TARFILE)

# Style variables
T= "	"
T2=$(T)$(T)
PREFIXNOTE="Must be an absolute path directory. Default: $(PWD)"

###################################

default: help

help:
	@echo 1. INSTALL CHOICES
	@echo $(T) all - Perform install, env, and clean. All parameters are valid to use here.
	@echo $(T) install - copy all files over to an installation directory. Installs most prerequisites.
	@echo $(T2) PREFIX=$(PREFIX) $(PREFIXNOTE)
	@echo $(T) cuttingedge - download and install the most up to date code. Does not include 'make env' or any prerequisites. Can be used instead of 'make install'
	@echo $(T2) PREFIX=$(PREFIX) $(PREFIXNOTE)
	@echo
	@echo 2. ENVIRONMENT CHECK
	@echo $(T) env - put all environmental variables into a profile file
	@echo $(T2) PROFILE=$(PROFILE)
	@echo $(T) check - check to see if all prerequisites are installed
	@echo $(T) test - create a test project using the test data found in the installation directory
	@echo
	@echo 3. OTHER
	@echo $(T) clean - delete the temporary files. Does not remove the result of 'make env.'
	@echo $(T2) PREFIX=$(PREFIX) $(PREFIXNOTE)
	@echo $(T2) NUMCPUS=$(NUMCPUS)
	@echo $(T2) PREFIX=$(PREFIX) $(PREFIXNOTE)
	@echo $(T2) PROJECT=$(PROJECT)

all: install env clean

install: install-prerequisites
	@echo "Don't forget to set up update PATH and PERL5LIB to $(PREFIX)/scripts and $(PREFIX)/lib"
	@echo "'make env' performs this step for you"
	@echo "DONE: installation of Lyve-SET v$(VERSION) complete."

install-prerequisites: install-mkdir install-vcftools install-CGP install-callsam install-SGELK install-varscan
	@echo DONE installing prerequisites

install-mkdir:
	-mkdir $(PREFIX)/build $(PREFIX)/lib $(PREFIX)/scripts

install-callsam:
	rm -rf $(PREFIX)/lib/callsam # make sure it's removed
	git clone https://github.com/lskatz/callsam.git $(PREFIX)/lib/callsam
install-SGELK:
	rm -rf $(TMPDIR)/Schedule # make sure it's removed
	git clone https://github.com/lskatz/Schedule--SGELK.git $(TMPDIR)/Schedule
	-mkdir -p $(PREFIX)/lib/Schedule
	mv -v $(TMPDIR)/Schedule/SGELK.pm $(PREFIX)/lib/Schedule/
	mv -v $(TMPDIR)/Schedule/README.md $(PREFIX)/lib/Schedule/
install-CGP:
	# make sure these scripts are removed
	rm -f $(PREFIX)/scripts/run_assembly_isFastqPE.pl
	rm -f $(PREFIX)/scripts/run_assembly_trimClean.pl
	rm -f $(PREFIX)/scripts/run_assembly_shuffleReads.pl
	rm -rvf $(PREFIX)/lib/cg-pipeline-code
	# CGP scripts that are needed and that don't depend on CGP libraries
	svn checkout https://svn.code.sf.net/p/cg-pipeline/code/ $(PREFIX)/lib/cg-pipeline-code
	ln -s $(PREFIX)/lib/cg-pipeline-code/cg_pipeline/branches/lkatz/scripts/run_assembly_isFastqPE.pl $(PREFIX)/scripts/
	ln -s $(PREFIX)/lib/cg-pipeline-code/cg_pipeline/branches/lkatz/scripts/run_assembly_trimClean.pl $(PREFIX)/scripts/
	ln -s $(PREFIX)/lib/cg-pipeline-code/cg_pipeline/branches/lkatz/scripts/run_assembly_shuffleReads.pl $(PREFIX)/scripts/

install-vcftools:
	rm -rvf $(PREFIX)/lib/vcftools_0.1.12b
	rm -vf $(PREFIX)/scripts/vcf-sort
	rm -vf $(PREFIX)/lib/Vcf.pm
	wget 'http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F&ts=1409260024&use_mirror=ufpr' -O $(TMPDIR)/vcftools_0.1.12b.tar.gz
	cd $(TMPDIR) && \
	tar zxvf vcftools_0.1.12b.tar.gz
	mv $(TMPDIR)/vcftools_0.1.12b $(PREFIX)/lib/
	cd $(PREFIX)/lib/vcftools_0.1.12b &&\
    make --directory=$(PREFIX)/lib/vcftools_0.1.12b MAKEFLAGS=""
	ln -s $(PREFIX)/lib/vcftools_0.1.12b/perl/vcf-sort $(PREFIX)/scripts/
	ln -s $(PREFIX)/lib/vcftools_0.1.12b/perl/Vcf.pm $(PREFIX)/lib/

install-varscan:
	wget 'http://downloads.sourceforge.net/project/varscan/VarScan.v2.3.7.jar?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvarscan%2Ffiles%2F&ts=1413398147&use_mirror=ufpr' -O $(PREFIX)/lib/varscan.v2.3.7.jar

cuttingedge: install-mkdir cuttingedge-gitclone install-prerequisites
	@echo "DONE installing the cutting edge version"

cuttingedge-gitclone:
	git clone --recursive https://github.com/lskatz/lyve-SET.git $(PREFIX)/build/Lyve-SET
	rm -r $(PREFIX)/build/Lyve-SET/build # avoid directory-not-empty error
	mv -vf $(PREFIX)/build/Lyve-SET/* $(PREFIX)/


env:
	echo "#Lyve-SET" >> $(PROFILE)
	echo "export PATH=\$$PATH:$(PREFIX)/scripts" >> $(PROFILE)
	echo "export PERL5LIB=\$$PERL5LIB:$(PREFIX)/lib" >> $(PROFILE)

clean:
	rm -vrf $(TMPDIR)/*
	@echo "Remember to remove the line with PATH and Lyve-SET from $(PROFILE)"

test:
	@echo "Test data set given by CFSAN's snp-pipeline package found at https://github.com/CFSAN-Biostatistics/snp-pipeline"
	set_manage.pl --create $(PROJECT)
	set_manage.pl $(PROJECT) --add-reads $(PREFIX)/testdata/reads/sample1.fastq.gz
	set_manage.pl $(PROJECT) --add-reads $(PREFIX)/testdata/reads/sample2.fastq.gz
	set_manage.pl $(PROJECT) --add-reads $(PREFIX)/testdata/reads/sample3.fastq.gz
	set_manage.pl $(PROJECT) --add-reads $(PREFIX)/testdata/reads/sample4.fastq.gz
	set_manage.pl $(PROJECT) --change-reference $(PREFIX)/testdata/reference/lambda_virus.fasta
	set_manage.pl $(PROJECT) --add-assembly $(PREFIX)/testdata/reference/lambda_virus.fasta
	launch_set.pl $(PROJECT) --noclean --snpcaller callsam --msa-creation lyve-set-lowmem --numcpus $(NUMCPUS)

check: check-sys check-Lyve-SET-PATH check-CGP-assembly check-Lyve-SET check-PERL check-smalt check-freebayes check-raxml check-freebayes check-phyml
	@echo --OK
check-sys:
	@F=$$(which cat) && echo "Found $$F"
	@F=$$(which gzip) && echo "Found $$F"
	@F=$$(which perl) && echo "Found $$F"
check-smalt:
	@F=$$(which smalt 2>/dev/null) && echo "Found smalt at $$F"
check-freebayes:
	@F=$$(which freebayes 2>/dev/null) && echo "Found freebayes at $$F"
check-raxml:
	@ RAXML=$$((which raxml || which raxmlHPC-PTHREADS) 2>/dev/null) && echo "Found raxml at $$RAXML"
check-phyml:
	@ PHYML=$$((which phyml || which phyml_linux_64 ) 2>/dev/null) && echo "Found phyml at $$PHYML"
check-callsam:
	@export PATH=$$PATH:$(PREFIX)/lib/callsam/bin && which callsam_MT.pl >/dev/null
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
	@echo "Looking for Graph::Centrality::Pagerank"
	@perl -I $(PREFIX)/lib -MGraph::Centrality::Pagerank -e 1

