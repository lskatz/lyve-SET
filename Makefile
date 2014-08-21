# Author: Lee Katz <lkatz@cdc.gov>
# Lyve-SET

PREFIX := /opt/Lyve-SET
PROFILE := $(HOME)/.bashrc
VERSION := 0.8.1
PROJECT := "setTestProject"

# Derived variables
TMPDIR := $(PREFIX)/build
TARFILE=Lyve-SET.v$(VERSION).tar.gz
TMPTARFILE=$(TMPDIR)/$(TARFILE)

# Style variables
T= "	"
T2=$(T)$(T)

###################################

default: help

help:
	@echo Commands:
	@echo $(T) all - Perform install, env, and clean. All parameters are valid to use here.
	@echo $(T) install - copy all files over to an installation directory
	@echo $(T2) PREFIX=$(PREFIX)
	@echo $(T2) VERSION=$(VERSION)
	@echo $(T) cuttingedge - download and install the most up to date code. Does not include 'make env'
	@echo $(T2) PREFIX=$(PREFIX)
	@echo $(T) env - put all environmental variables into a profile file 
	@echo $(T2) PROFILE=$(PROFILE)
	@echo $(T) clean - delete the temporary files. Does not remove the result of 'make env.'
	@echo $(T2) PREFIX=$(PREFIX)
	@echo $(T) test - create a test project using the test data found in the installation directory
	@echo $(T2) PREFIX=$(PREFIX)
	@echo $(T2) PROJECT=$(PROJECT)
	@echo NOTES: 
	@echo $(T) All paths must be absolute
	@echo Example:
	@echo $(T) make all PREFIX=$(PREFIX) VERSION=$(VERSION) PROFILE=$(PROFILE)
	@echo $(T) "make cuttingedge PREFIX=$(PREFIX) && make env PROFILE=$(PROFILE) && make test PREFIX=$(PREFIX) PROJECT=$(PROJECT)"

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
	cd $(PREFIX)/lib && \
	git clone https://github.com/lskatz/callsam.git
	git clone https://github.com/lskatz/Schedule--SGELK.git $(PREFIX)/lib/Schedule

cuttingedge:
	git clone --recursive https://github.com/lskatz/lyve-SET.git $(PREFIX)

env:
	echo -e "#Lyve-SET\nexport PATH=\$$PATH:$(PREFIX)" >> $(PROFILE)

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
	launch_set.pl $(PROJECT) --noclean --snpcaller callsam --msa-creation lyve-set-lowmem

fail:
	touch /dfjkd/dfjdksajo/dfj32098/dkdl
	exit 5
	exit 1
