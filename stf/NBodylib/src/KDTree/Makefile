include ../../../Makefile.config

all : libKD.a

CFILES = $(wildcard *.cxx)
INCL = $(wildcard *.h)
OBJ = $(patsubst %.cxx,%.o,$(wildcard *.cxx))
MAKECHECKS = Makefile ../../../Makefile.config
LIBCHECKS = ../Math/libMath.a ../NBody/libNBody.a

LIB=libKD.a

$(LIB) :  $(INCL) $(OBJ)
	ar r $(LIB) *.o 
	cp $(LIB) $(NBODYLIBDIR)
	cp *.h $(NBODYINCLUDEDIR)

%.o : %.cxx $(INCL) $(MAKECHECKS) $(LIBCHECKS)
	$(C+) $(NBODYC+FLAGS) $(NBODYIFLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.a
