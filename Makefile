CC     = cc -ansi -pedantic -Wall
COPT   = -O3 -L$(HOME)/lib -I$(HOME)/include
OFILES = genloop.o
LIBS   = -lbiop -lgen -lm -lxml2

genloop : $(OFILES)
	$(CC) $(COPT) -o $@ $(OFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f $(OFILES)
