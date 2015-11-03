CFILES=$(wildcard src/*/*.cxx) 
HFILES = $(wildcard src/*/*.h)

all:
	cd halotree; make 
	cd baryons; make 

clean:
	cd halotree; make clean
	cd baryons; make clean


ifeq "$(wildcard doc/doxy.log)" ""
doc: doc_
doc_:
	cd doc; doxygen Doxyfile > doxy.log;
else
doc: doc/Doxyfile $(CFILES) $(HFILES)
	cd doc; doxygen Doxyfile > doxy.log;
endif

docclean: docc_
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
