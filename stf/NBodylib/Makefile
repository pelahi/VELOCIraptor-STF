CFILES=$(wildcard src/*/*.cxx)
HFILES = $(wildcard src/*/*.h)

all:
	mkdir -p lib; mkdir -p include; cd src; make

clean:
	cd src; make clean;


ifeq "$(wildcard doc/doxy.log)" ""
doc: doc_
doc_:
	cd doc; mkdir -p html; mkdir -p latex; mkdir -p xml; doxygen Doxyfile > doxy.log;
else
doc: doc/Doxyfile $(CFILES) $(HFILES)
	cd doc; mkdir -p html; mkdir -p latex; mkdir -p xml; doxygen Doxyfile > doxy.log;
endif

docclean: docc_
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
