CC = cc
COPT = -ansi -pedantic -Wall -L$(HOME)/lib -I$(HOME)/include -O3
LIB  = -lbiop -lgen -lm -lxml2

EXE  = topscan mergestride

all : $(EXE)

topscan : topscan.o
	$(CC) $(COPT) -o $@ $< $(LIB)

mergestride : mergestride.o
	$(CC) $(COPT) -o $@ $< $(LIB)

.c.o :
	$(CC) $(COPT) -c -o $@ $<
