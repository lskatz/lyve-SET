# Author: Lee Katz <lkatz@cdc.gov>
# Lyve-SET
 
PREFIX := $(PWD)
PROFILE := $(HOME)/.bashrc
PROJECT := "setTestProject"
NUMCPUS := 1
SHELL   := /bin/bash

# Derived variables
TMPDIR := $(PREFIX)/build

# Style variables
T= "	"
T2=$(T)$(T)
PREFIXNOTE="Must be an absolute path directory. Default: $(PWD)"

###################################

default: install 

help:
	@echo "Please see README.md for additional help"

all: install env clean

install: install-prerequisites
	@echo "Don't forget to include $(PREFIX)/scripts in your PATH"
	@echo "DONE: installation of Lyve-SET complete."

install-prerequisites: lib/vcftools_0.1.12b lib/cg-pipeline lib/Schedule/SGELK.pm lib/varscan.v2.3.7.jar lib/phast/phast.faa lib/samtools-1.2/samtools lib/bcftools-1.2/bcftools lib/smalt-0.7.6/bin/smalt lib/snap/snap-aligner lib/standard-RAxML-8.1.16/raxmlHPC install-perlModules install-config lib/snpEff.jar lib/stampy/stampy.py
	@echo DONE installing prerequisites

clean: clean-vcftools clean-CGP clean-SGELK clean-varscan clean-phast clean-samtools clean-bcftools clean-smalt clean-snap clean-raxml clean-perlModules clean-config clean-snpEff clean-stampy
	@echo "Remember to remove the Lyve-SET path from PATH"

lib/vcftools_0.1.12b:
	wget 'http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F&ts=1409260024&use_mirror=ufpr' -O $(TMPDIR)/vcftools_0.1.12b.tar.gz
	cd $(TMPDIR) && \
	tar zxvf vcftools_0.1.12b.tar.gz
	mv $(TMPDIR)/vcftools_0.1.12b $(PREFIX)/lib/
	cd $(PREFIX)/lib/vcftools_0.1.12b &&\
    make --directory=$(PREFIX)/lib/vcftools_0.1.12b MAKEFLAGS=""
	ln -s $(PREFIX)/lib/vcftools_0.1.12b/perl/Vcf.pm $(PREFIX)/lib/
	ln -s $(PREFIX)/lib/vcftools_0.1.12b/perl/vcf-sort $(PREFIX)/scripts/

clean-vcftools:
	rm -rvf $(PREFIX)/lib/vcftools_0.1.12b
	rm -vf $(PREFIX)/scripts/vcf-sort
	rm -vf $(PREFIX)/lib/Vcf.pm

lib/cg-pipeline:
	# CGP scripts that are needed and that don't depend on CGP libraries
	git clone https://github.com/lskatz/cg-pipeline $(PREFIX)/lib/cg-pipeline
	ln -s $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_isFastqPE.pl $(PREFIX)/scripts/
	ln -s $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_trimClean.pl $(PREFIX)/scripts/
	ln -s $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_shuffleReads.pl $(PREFIX)/scripts/
	ln -s $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_removeDuplicateReads.pl $(PREFIX)/scripts/
	ln -s $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_readMetrics.pl $(PREFIX)/scripts/
	ln -s $(PREFIX)/lib/cg-pipeline/scripts/run_assembly_metrics.pl $(PREFIX)/scripts/

clean-CGP:
	rm -rvf $(PREFIX)/lib/cg-pipeline
	rm -vf $(PREFIX)/scripts/run_assembly_*

lib/Schedule/SGELK.pm:
	git clone https://github.com/lskatz/Schedule--SGELK.git $(TMPDIR)/Schedule
	-mkdir -p $(PREFIX)/lib/Schedule
	mv -v $(TMPDIR)/Schedule/SGELK.pm $(PREFIX)/lib/Schedule/
	mv -v $(TMPDIR)/Schedule/README.md $(PREFIX)/lib/Schedule/
	mv -v $(TMPDIR)/Schedule/.git $(PREFIX)/lib/Schedule/

clean-SGELK:
	rm -rfv $(PREFIX)/lib/Schedule

lib/varscan.v2.3.7.jar:
	wget 'http://downloads.sourceforge.net/project/varscan/VarScan.v2.3.7.jar' -O $(PREFIX)/lib/varscan.v2.3.7.jar

clean-varscan:
	rm -vf $(PREFIX)/lib/varscan.v2.3.7.jar

lib/phast/phast.faa: 
	mkdir -p $(PREFIX)/lib/phast
	wget http://phast.wishartlab.com/phage_finder/DB/prophage_virus.db -O $(PREFIX)/$@
	makeblastdb -in $(PREFIX)/lib/phast/phast.faa -dbtype prot

clean-phast:
	rm -rvf $(PREFIX)/lib/phast

install-phispy:
	mkdir -p $(PREFIX)/lib/phispy
	wget http://downloads.sourceforge.net/project/phispy/phiSpyNov11_v2.3.zip -O $(PREFIX)/lib/phispy/phiSpyNov11_v2.3.zip
	cd $(PREFIX)/lib/phispy && unzip -o phiSpyNov11_v2.3.zip
	rm $(PREFIX)/lib/phispy/phiSpyNov11_v2.3.zip
	cd $(PREFIX)/lib/phispy/phiSpyNov11_v2.3 && make

clean-phispy:
	rm -rvf $(PREFIX)/lib/phispy

lib/samtools-1.2/samtools:
	wget 'https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2' -O $(TMPDIR)/samtools-1.2.tar.bz2
	cd $(TMPDIR) && tar jxvf samtools-1.2.tar.bz2
	mv $(TMPDIR)/samtools-1.2 $(PREFIX)/lib
	cd $(PREFIX)/lib/samtools-1.2 && make DFLAGS="-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_CURSES_LIB=0" LIBCURSES="" 
	cd $(PREFIX)/lib/samtools-1.2/htslib-1.2.1 && make
	ln -sf $(PREFIX)/lib/samtools-1.2/samtools $(PREFIX)/scripts/
	ln -sf $(PREFIX)/lib/samtools-1.2/misc/wgsim $(PREFIX)/scripts/
	ln -sf $(PREFIX)/lib/samtools-1.2/htslib-1.2.1/bgzip $(PREFIX)/scripts
	ln -sf $(PREFIX)/lib/samtools-1.2/htslib-1.2.1/tabix $(PREFIX)/scripts

clean-samtools:
	rm -rvf $(PREFIX)/lib/samtools*
	rm -f $(PREFIX)/scripts/bgzip $(PREFIX)/scripts/tabix $(PREFIX)/scripts/wgsim $(PREFIX)/scripts/samtools

lib/bcftools-1.2/bcftools:
	wget 'https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2' -O $(TMPDIR)/bcftools-1.2.tar.bz2
	cd $(TMPDIR) && tar jxvf bcftools-1.2.tar.bz2
	mv $(TMPDIR)/bcftools-1.2 $(PREFIX)/lib/bcftools-1.2
	cd $(PREFIX)/lib/bcftools-1.2 && make
	ln -sf $(PREFIX)/lib/bcftools-1.2/bcftools $(PREFIX)/scripts
	ln -sf $(PREFIX)/lib/bcftools-1.2/vcfutils.pl $(PREFIX)/scripts

clean-bcftools:
	rm -rfv $(PREFIX)/lib/bcftools*
	rm -f   $(PREFIX)/lib/bcftools-1.2/bcftools $(PREFIX)/lib/bcftools-1.2/vcfutils.pl

lib/smalt-0.7.6/bin/smalt:
	wget --continue 'http://downloads.sourceforge.net/project/smalt/smalt-0.7.6-static.tar.gz' -O $(TMPDIR)/smalt-0.7.6-static.tar.gz
	cd $(TMPDIR) && tar zxvf smalt-0.7.6-static.tar.gz
	mv $(TMPDIR)/smalt-0.7.6 $(PREFIX)/lib/
	cd $(PREFIX)/lib/smalt-0.7.6 && ./configure --prefix $(PREFIX)/lib/smalt-0.7.6 && make && make install
	for i in basqcol fetchseq mixreads readstats simqual simread smalt splitmates splitreads trunkreads; do \
	  ln -svf $(PREFIX)/lib/smalt-0.7.6/bin/$$i $(PREFIX)/scripts/; \
	done

clean-smalt:
	rm -rvf $(PREFIX)/lib/smalt*
	for i in basqcol fetchseq mixreads readstats simqual simread smalt splitmates splitreads trunkreads; do \
	  rm -fv $(PREFIX)/scripts/$$i; \
	done

lib/snap/snap-aligner:
	git clone https://github.com/amplab/snap.git $(PREFIX)/lib/snap
	cd $(PREFIX)/lib/snap && git checkout v1.0beta.18 && make
	cp -vn $(PREFIX)/lib/snap/snap-aligner  $(PREFIX)/scripts/snap
	cd $(PREFIX)/lib/snap && make clean && make CXXFLAGS="-DLONG_READS -O3 -Wno-format -MMD -ISNAPLib -msse"
	cp -vn $(PREFIX)/lib/snap/snap-aligner  $(PREFIX)/scripts/snapxl

clean-snap:
	rm -rvf $(PREFIX)/lib/snap 
	rm -f   $(PREFIX)/scripts/snap $(PREFIX)/scripts/snapxl

lib/standard-RAxML-8.1.16/raxmlHPC:
	wget 'https://github.com/stamatak/standard-RAxML/archive/v8.1.16.tar.gz' -O $(TMPDIR)/raxml_v8.1.16.tar.gz
	cd $(TMPDIR) && tar zxvf raxml_v8.1.16.tar.gz
	mv $(TMPDIR)/standard-RAxML-8.1.16 $(PREFIX)/lib
	cd $(PREFIX)/lib/standard-RAxML-8.1.16 && (for i in Makefile.*; do make -f $$i; rm -vf *.o; done;)
	find $(PREFIX)/lib/standard-RAxML-8.1.16 -type f -executable -exec ln -sv {} $(PREFIX)/scripts \;

clean-raxml:
	rm -rvf $(PREFIX)/lib/standard-RAxML-8.1.16
	rm -fv  $(PREFIX)/scripts/raxmlHPC*

# There isn't an easy way for Make to detect whether these files have
# been installed and so it shouldn't depend on a file existing and
# should run no matter what. CPAN will know if the module will be
# there.
# TODO I should probably have an individual install command per 
# module so that Make can do its job.
install-perlModules:
	@echo "Installing Perl modules using cpanminus"
	for package in Config::Simple File::Slurp Math::Round Number::Range Statistics::Distributions Statistics::Basic Graph::Centrality::Pagerank String::Escape Statistics::LineFit; do \
	  perl scripts/cpanm --self-contained -L $(PREFIX)/lib $$package; \
		if [ $$? -gt 0 ]; then exit 1; fi; \
	done;
	@echo "Done with Perl modules"

clean-perlModules:
	@echo "Perl modules were installed using CPAN which doesn't have an uninstalling mechanism"

install-config: config/LyveSET.conf config/presets.conf

config/%.conf:
	b=`basename $@` && d=`dirname $@` &&\
	  cp -v config/original/$$b config/$$b

clean-config:
	rm -v $(PREFIX)/config/*.conf

lib/snpEff.jar:
	wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip -O $(PREFIX)/lib/snpEff_latest_core.zip
	cd $(PREFIX)/lib && unzip -o snpEff_latest_core.zip
	cp $(PREFIX)/lib/snpEff/snpEff.jar $(PREFIX)/lib/
	cp $(PREFIX)/lib/snpEff/snpEff.config $(PREFIX)/config/snpEff.conf
	rm -rf $(PREFIX)/lib/snpEff_latest_core.zip
	rm -rf $(PREFIX)/lib/snpEff

clean-snpEff:
	rm -rvf $(PREFIX)/lib/snpEff.jar

lib/stampy/stampy.py:
	wget www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz -O $(TMPDIR)/Stampy-latest.tgz
	mkdir $(PREFIX)/lib/stampy && tar zxvf $(TMPDIR)/Stampy-latest.tgz -C $(PREFIX)/lib/stampy --strip-components 1
	cd $(PREFIX)/lib/stampy && make CXXFLAGS='-Wno-deprecated -Wunused-but-set-variable'
	ln -sv $(PREFIX)/lib/stampy/stampy.py $(PREFIX)/scripts

clean-stampy:
	rm -rvf $(PREFIX)/lib/stampy
	rm $(PREFIX)/scripts/stampy.py

env:
	echo "#Lyve-SET" >> $(PROFILE)
	echo "export PATH=\$$PATH:$(PREFIX)/scripts" >> $(PROFILE)

check: check-sys check-Lyve-SET-PATH check-CGP-assembly check-Lyve-SET check-PERL check-smalt check-freebayes check-raxml check-freebayes check-blast
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
check-blast:
	which blastx
	which blastp
	which makeblastdb
