include ../../../Makefile.config

all : libInitCond.a

CFILES = $(wildcard *.cxx)
INCL = $(wildcard *.h)
OBJ = $(patsubst %.cxx,%.o,$(wildcard *.cxx))
MAKECHECKS = Makefile ../../../Makefile.config
LIBCHECKS = ../Math/libMath.a ../NBody/libNBody.a

libInitCond.a : $(INCL) $(OBJ)
	ar r libInitCond.a *.o
	cp libInitCond.a $(NBODYLIBDIR)
	cp *.h $(NBODYINCLUDEDIR)

%.o : %.cxx %.h $(MAKECHECKS) $(LIBCHECKS)
	$(C+) $(NBODYC+FLAGS) $(NBODYIFLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.a
