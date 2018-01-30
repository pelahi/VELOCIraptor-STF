include ../../Makefile.config
MAKECHECK=../../Makefile.config Makefile

OBJS = $(patsubst %.cxx,%.o,$(wildcard *.cxx))
INCL   = *.h
EXEC = treefrog

all : $(EXEC)

$(EXEC) : $(OBJS)
	$(C+) -o $(EXEC) $(C+FLAGS) $(OBJS) $(LFLAGS) $(C+LIBS)
	cp $(EXEC) $(STFBINDIR)

%.o: %.cxx $(INCL) $(MAKECHECK) $(LIBCHECK)
	$(C+) $(C+FLAGS) $(IFLAGS) -c -o $@ $<

.PHONY : clean

clean :
	rm -f $(OBJS) $(EXEC)

ifeq "$(wildcard doc/doxy.log)" ""
doc: doc_
doc_:
	cd doc; mkdir -p html; mkdir -p latex; mkdir -p xml; doxygen Doxyfile > doxy.log;
else
doc:  doc/Doxyfile $(CFILES) $(HFILES)
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
