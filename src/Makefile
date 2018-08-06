DEST   = $(HOME)/bin
LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include
CC     = cc
COPT   = -ansi -pedantic -Wall -L$(LIBDIR) -I$(INCDIR) -O3
LIB    = -lbiop -lgen -lm -lxml2
EXE    = topscan mergestride

all : $(EXE)

topscan : topscan.o
	$(CC) $(COPT) -o $@ $< $(LIB)

mergestride : mergestride.o
	$(CC) $(COPT) -o $@ $< $(LIB)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f topscan.o mergestride.o

distclean : clean
	\rm -f $(EXE)

install :
	cp $(EXE) $(DEST)
