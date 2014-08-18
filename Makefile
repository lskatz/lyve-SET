# Author: Lee Katz <lkatz@cdc.gov>
# Lyve-SET

PREFIX := /opt/Lyve-SET
PROFILE := $(HOME)/.bashrc
VERSION := 0.8

TMPDIR := $(PREFIX)/build
TARFILE=Lyve-SET.v$(VERSION).tar.gz
TMPTARFILE=$(TMPDIR)/$(TARFILE)

T= "	"
T2=$(T)$(T)

default: help

help:
	@echo Commands:
	@echo $(T) all - Perform install and env. All parameters are valid to use here.
	@echo $(T) install - copy all files over to an installation directory
	@echo $(T2) PREFIX=$(PREFIX) 
	@echo $(T2) VERSION=$(VERSION)
	@echo $(T) env - put all environmental variables into a profile file 
	@echo $(T2) PROFILE=$(PROFILE)
	@echo Example:
	@echo $(T) make all PREFIX=$(PREFIX) VERSION=$(VERSION) PROFILE=$(PROFILE)

all: install env

install:
	mkdir $(PREFIX) 
	mkdir $(TMPDIR)
	wget https://github.com/lskatz/lyve-SET/archive/v$(VERSION).tar.gz -O $(TMPTARFILE)
	cd $(TMPDIR) && \
	tar zxvf $(TMPTARFILE)
	# Move all the untarred files to the install directory
	mv -v $(TMPDIR)/lyve-SET-$(VERSION)/* $(PREFIX)/
	@echo NOTE: 'make env' in order to set your path permanently.
	@echo NOTE: 'make clean' to remove the temporary directory.

env:
	echo -e "#Lyve-SET\nexport PATH=\$$PATH:$(PREFIX)" >> $(PROFILE)

clean:
	rm -vr $(TMPDIR)
	@echo "remove the line with PATH and Lyve-SET from $(PROFILE)"

