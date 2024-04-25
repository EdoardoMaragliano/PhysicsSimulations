CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall 
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`



cristallo.x: main.cpp Reticolo.o
	$(CC) $(CFLAGS) -o main.x main.cpp Reticolo.o $(CFLAGSROOT) $(LIBSROOT)
Reticolo.o: Reticolo.hpp Reticolo.cpp
	$(CC) $(CFLAGS) -c Reticolo.cpp $(CFLAGSROOT) $(LIBSROOT) 

clean:
	rm *.o 