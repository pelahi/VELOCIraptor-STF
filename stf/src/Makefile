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
	ar -x $(NBODYLIBDIR)/libNBody.a
	ar r $(LIB) *.o
	cp $(LIB) $(STFLIBDIR)
	cp $(INCL) $(STFINCLUDEDIR)

%.o: %.cxx $(INCL) $(MAKECHECK) $(LIBCHECK) 
	$(C+) $(C+FLAGS) $(IFLAGS) -c -o $@ $< 

.PHONY : clean

clean :
	rm -f $(OBJS) $(EXEC)
	rm -f *.a
