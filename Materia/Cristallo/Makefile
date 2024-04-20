CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall 
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`



cristallo.x: cristallo.cpp Reticolo.o
	$(CC) $(CFLAGS) -o cristallo.x cristallo.cpp Reticolo.o $(CFLAGSROOT) $(LIBSROOT)
Reticolo.o: Reticolo.hpp Reticolo.cpp
	$(CC) $(CFLAGS) -c Reticolo.cpp $(CFLAGSROOT) $(LIBSROOT) 

clean:
	rm *.o 