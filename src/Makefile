LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include
CC     = cc
COPT   = -ansi -pedantic -Wall -L$(LIBDIR) -I$(INCDIR) -g
LIB    = -lbiop -lgen -lm -lxml2
EXE    = topscan mergestride mergepdbsecstr

all : $(EXE)

topscan : topscan.o
	$(CC) $(COPT) -o $@ $< $(LIB)

mergestride : mergestride.o
	$(CC) $(COPT) -o $@ $< $(LIB)

mergepdbsecstr : mergepdbsecstr.o
	$(CC) $(COPT) -o $@ $< $(LIB)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f topscan.o mergestride.o mergepdbsecstr.o

distclean : clean
	\rm -f $(EXE)

install :
	./install.sh

test :
	(cd t; ./runtest.sh)

