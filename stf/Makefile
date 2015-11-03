include Makefile.config
all: NBodylib_ stf_ analysis_
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

CFILES = $(wildcard src/*.cxx) 
HFILES = $(wildcard src/*.h)

binc_ : 
	cd bin; rm *
NBodylib_ : Makefile.config
	cd NBodylib; make
NBodylibc_ : 
	cd NBodylib; make clean
stf_ : Makefile.config
	cd src; make
stfc_ : 
	cd src; make clean

analysis_ : Makefile.config
	cd analysis; make
analysisc_ : 
	cd analysis; make clean

ifeq "$(wildcard doc/doxy.log)" ""
doc: doc_
doc_:
	cd NBodylib; make doc;
	cd doc; doxygen Doxyfile > doxy.log;
else
doc:  doc/Doxyfile $(CFILES) $(HFILES)
	cd NBodylib; make doc;
	cd doc; doxygen Doxyfile > doxy.log;
endif
docc_:
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

