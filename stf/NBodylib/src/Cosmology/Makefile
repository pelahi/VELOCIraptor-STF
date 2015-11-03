include ../../../Makefile.config

all : libCosmology.a

CFILES = $(wildcard *.cxx)
INCL = $(wildcard *.h)
OBJ = $(patsubst %.cxx,%.o,$(wildcard *.cxx))
MAKECHECKS = Makefile ../../../Makefile.config
LIBCHECKS = ../Math/libMath.a

libCosmology.a : $(INCL) $(OBJ)
	ar r libCosmology.a Cosmology.o
	cp libCosmology.a $(NBODYLIBDIR)
	cp *.h $(NBODYINCLUDEDIR)

%.o : %.cxx %.h $(MAKECHECKS) $(LIBCHECKS)
	$(C+) $(NBODYC+FLAGS) $(NBODYIFLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.a
