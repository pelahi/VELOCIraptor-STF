include ../../../Makefile.config

all : libMath.a

CFILES = $(wildcard *.cxx)
INCL = $(wildcard *.h)
OBJ = $(patsubst %.cxx,%.o,$(wildcard *.cxx))
MAKECHECKS = Makefile ../../../Makefile.config
LIBCHECKS = 

LIB = libMath.a

$(LIB) : $(INCL) $(OBJ)
	ar r $(LIB) *.o 
	cp $(LIB) $(NBODYLIBDIR)
	cp *.h $(NBODYINCLUDEDIR)

%.o : %.cxx %.h $(MAKECHECKS) $(LIBCHECKS)
	$(C+) $(NBODYC+FLAGS) $(NBODYIFLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.a
