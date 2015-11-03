include ../Makefile.config
MAKECHECK=../Makefile.config Makefile

OBJS = $(patsubst %.cxx,%.o,$(wildcard *.cxx))
INCL   = *.h

all : $(EXEC)

$(EXEC) : $(OBJS)
	$(C+) -o $(EXEC) $(C+FLAGS) $(OBJS) $(LFLAGS) $(C+LIBS)  
	cp $(EXEC) $(STFBINDIR)

%.o: %.cxx $(INCL) $(MAKECHECK) $(LIBCHECK) 
	$(C+) $(C+FLAGS) $(IFLAGS) -c -o $@ $< 

.PHONY : clean

clean :
	rm -f $(OBJS) $(EXEC)
