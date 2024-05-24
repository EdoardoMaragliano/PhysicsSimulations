CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall 
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`



cristallo.x: main.cpp Lattice.o
	$(CC) $(CFLAGS) -o main.x main.cpp Lattice.o $(CFLAGSROOT) $(LIBSROOT)
Lattice.o: Lattice.hpp Lattice.cpp
	$(CC) $(CFLAGS) -c Lattice.cpp $(CFLAGSROOT) $(LIBSROOT) 

clean:
	rm *.o 