# Author: Lee Katz <lkatz@cdc.gov>
# Lyve-SET

PREFIX := /opt/Lyve-SET
PROFILE := $(HOME)/.bashrc
VERSION := 0.8.3
PROJECT := "setTestProject"
NUMCPUS := 1
SHELL   := /bin/bash

# Derived variables
TMPDIR := $(PREFIX)/build
TARFILE=Lyve-SET.v$(VERSION).tar.gz
TMPTARFILE=$(TMPDIR)/$(TARFILE)

# Style variables
T= "	"
T2=$(T)$(T)
PREFIXNOTE="Must be an absolute path directory"

###################################

default: help

help:
	@echo 1. INSTALL CHOICES
	@echo $(T) all - Perform install, env, and clean. All parameters are valid to use here.
	@echo $(T) install - copy all files over to an installation directory
	@echo $(T2) PREFIX=$(PREFIX) $(PREFIXNOTE)
	@echo $(T2) VERSION=$(VERSION)
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

install:
	mkdir $(PREFIX) 
	mkdir $(TMPDIR)
	wget https://github.com/lskatz/lyve-SET/archive/v$(VERSION).tar.gz -O $(TMPTARFILE)
	cd $(TMPDIR) && \
	tar zxvf $(TARFILE)
	# Move all the untarred files to the install directory
	mv -v $(TMPDIR)/lyve-SET-$(VERSION)/* $(PREFIX)/
	# download necessary submodules because git doesn't package them in the release
	rm -rvf $(PREFIX)/lib/*
	# Git submodules
	git clone https://github.com/lskatz/callsam.git $(PREFIX)/lib/callsam
	git clone https://github.com/lskatz/Schedule--SGELK.git $(PREFIX)/lib/Schedule
	# CGP scripts that are needed and that don't depend on CGP libraries
	svn checkout https://svn.code.sf.net/p/cg-pipeline/code/ $(PREFIX)/lib/cg-pipeline-code
	ln -s $(PREFIX)/lib/cg-pipeline-code/cg_pipeline/branches/lkatz/scripts/run_assembly_isFastqPE.pl $(PREFIX)/
	ln -s $(PREFIX)/lib/cg-pipeline-code/cg_pipeline/branches/lkatz/scripts/run_assembly_trimClean.pl $(PREFIX)/
	ln -s $(PREFIX)/lib/cg-pipeline-code/cg_pipeline/branches/lkatz/scripts/run_assembly_shuffleReads.pl $(PREFIX)/
	# vcftools
	wget 'http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F&ts=1409260024&use_mirror=ufpr' -O $(TMPDIR)/vcftools_0.1.12b.tar.gz
	cd $(TMPDIR) && \
	tar zxvf vcftools_0.1.12b.tar.gz
	mv $(TMPDIR)/vcftools_0.1.12b $(PREFIX)/lib/
	cd $(PREFIX)/lib/vcftools_0.1.12b && make
	ln -s $(PREFIX)/lib/vcftools_0.1.12b/bin/vcf-sort $(PREFIX)/
	ln -s $(PREFIX)/lib/vcftools_0.1.12b/lib/Vcf.pm $(PREFIX)/lib/

cuttingedge:
	git clone --recursive https://github.com/lskatz/lyve-SET.git $(PREFIX)

env:
	echo "#Lyve-SET" >> $(PROFILE)
	echo "export PATH=\$$PATH:$(PREFIX)" >> $(PROFILE)
	echo "export PERL5LIB=\$$PERL5LIB:$(PREFIX/lib)" >> $(PROFILE)

clean:
	rm -vrf $(TMPDIR)
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

check: check-cat check-gzip check-CGP-assembly check-Lyve-SET check-PERL
	@echo --OK
check-cat:
	@which cat >/dev/null
check-gzip:
	@which gzip >/dev/null
check-smalt:
	@which smalt >/dev/null
check-CGP-assembly:
	@which run_assembly_shuffleReads.pl run_assembly_trimClean.pl run_assembly_isFastqPE.pl >/dev/null
check-Lyve-SET:
	@echo Checking that SET is in your path
	@which launch_set.pl >/dev/null
check-PERL:
	@echo Checking for perl multithreading
	@perl -Mthreads -e 1
	@echo Checking for perl modules
	@perl -MFile::Slurp -e 1
	@perl -MString::Escape -e 1
	@perl -MGraph::Centrality::Pagerank -e 1
	
fail:
	touch /dfjkd/dfjdksajo/dfj32098/dkdl
	exit 5
	exit 1
