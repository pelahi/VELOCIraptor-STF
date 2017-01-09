CFILES=$(wildcard src/*/*.cxx) 
HFILES = $(wildcard src/*/*.h)

all:
	cd treefrog; make 
	cd baryons; make 

clean:
	cd treefrog; make clean
	cd baryons; make clean


ifeq "$(wildcard doc/doxy.log)" ""
doc: doc_
doc_:
	cd treefrog; make doc;
else
doc: doc/Doxyfile $(CFILES) $(HFILES)
	cd treefrog; make doc;
endif

docclean: docc_
docc_:
	cd treefrog; make docclean;
