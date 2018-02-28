include Makefile.config
all:
	make NBodylib_
	make stf_
	make analysis_
NBody: NBodylib_
lib: NBodylib_
analysis: analysis_
stf: stf_
clean: stfc_ binc_ analysisc_
libclean: NBodylibc_
analysisclean: analysisc_
binclean: binc_
docclean: docc_
allclean: stfc_ binc_ analysisc_ NBodylibc_ docc_
libstf: NBodylib_ libstf_

CFILES = $(wildcard src/*.cxx)
HFILES = $(wildcard src/*.h)
EXEFILES = $(wildcard bin/*)

binc_ :
	if [ -d "bin" ]; then for file in $(EXEFILES); do rm $$file; done; fi
NBodylib_ : Makefile.config
	cd NBodylib; make
NBodylibc_ :
	cd NBodylib; make clean
stf_ : Makefile.config
	mkdir -p bin
	cd src; make
stfc_ :
	cd src; make clean

analysis_ : Makefile.config
	mkdir -p bin
	cd analysis; make
analysisc_ :
	cd analysis; make clean

libstf_:
	mkdir -p $(STFLIBDIR)
	mkdir -p $(STFINCLUDEDIR)
	cd src; make libstf

ifeq "$(wildcard doc/doxy.log)" ""
doc: doc_
doc_:
	cd NBodylib; make doc;
	cd analysis; make doc;
	cd doc; mkdir -p html; mkdir -p latex; mkdir -p xml; doxygen Doxyfile > doxy.log;
else
doc:  doc/Doxyfile $(CFILES) $(HFILES)
	cd NBodylib; make doc;
	cd analysis; make doc;
	cd doc; mkdir -p html; mkdir -p latex; mkdir -p xml; doxygen Doxyfile > doxy.log;
endif
docc_:
ifneq "$(wildcard doc/xml)" ""
	rm -r doc/xml/
endif
ifneq "$(wildcard doc/html)" ""
	rm -r doc/html/
endif
ifneq "$(wildcard doc/latex)" ""
	rm -r doc/latex/
endif
ifneq "$(wildcard doc/doxy.log)" ""
	rm -r doc/doxy.log
endif
	cd NBodylib; make docclean;
	cd analysis; make docclean;
