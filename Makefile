CC=g++
CFLAGS=-O9 -DNDEBUG

all: index

Basicrmq.o: Basicrmq.cpp
	$(CC) $(CFLAGS) -c Basicrmq.cpp

DFUDSrmq.o: DFUDSrmq.cpp
	$(CC) $(CFLAGS) -c DFUDSrmq.cpp

index: Basicrmq.o DFUDSrmq.o
	ar rc rmqFischer.a Basicrmq.o DFUDSrmq.o

test: 
	@$(CC) $(CFLAGS) rmqFischerDFUDS.cpp -o rmqFischerTest rmqFischer.a 

clean:
	-rm *~ *.o *.bak 
cleanall:
	-rm *~ *.o *.bak .a
