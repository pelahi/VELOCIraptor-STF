include ../Makefile.config
MAKECHECK=../Makefile.config Makefile

OBJS = $(patsubst %.cxx,%.o,$(wildcard *.cxx))
OBJS_LIB = $(filter-out main.o,$(OBJS))
INCL = *.h
LIB = libstf.a

all : $(EXEC)
libstf : $(LIB)

$(EXEC) : $(OBJS)
	$(C+) -o $(EXEC) $(C+FLAGS) $(OBJS) $(LFLAGS) $(C+LIBS)
	cp $(EXEC) $(STFBINDIR)

$(LIB): $(INCL) $(OBJS_LIB)
	ar -x $(NBODYLIBDIR)/libAnalysis.a
	ar -x $(NBODYLIBDIR)/libCosmology.a
	ar -x $(NBODYLIBDIR)/libInitCond.a
	ar -x $(NBODYLIBDIR)/libKD.a
	ar -x $(NBODYLIBDIR)/libMath.a
	ar -x $(NBODYLIBDIR)/libNBody.a
	#ar -x $(GSL_DIR)/lib/libgsl.a
	#ar -x $(HDF_DIR)/lib/libhdf5.a
	#ar -x $(HDF_DIR)/lib/libhdf5_cpp.a
	ar rcs $(LIB) *.o
	cp $(LIB) $(STFLIBDIR)
	cp $(INCL) $(STFINCLUDEDIR)
	rm *.o

%.o: %.cxx $(INCL) $(MAKECHECK) $(LIBCHECK)
	$(C+) $(C+FLAGS) $(IFLAGS) -c -o $@ $<

.PHONY : clean

clean :
	rm -f $(OBJS) $(EXEC)
	rm -f *.a
