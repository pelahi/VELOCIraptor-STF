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
	cd halotree; make doc;
else
doc: doc/Doxyfile $(CFILES) $(HFILES)
	cd halotree; make doc;
endif

docclean: docc_
docc_:
	cd halotree; make docclean;
