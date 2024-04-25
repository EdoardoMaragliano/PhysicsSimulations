CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall 
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`



frattale.x: main.cpp Reticolo.o
	$(CC) $(CFLAGS) -o main.x main.cpp Reticolo.o $(CFLAGSROOT) $(LIBSROOT)
ReticoloF.o: ReticoloF.hpp Reticolo.cpp
	$(CC) $(CFLAGS) -c Reticolo.cpp $(CFLAGSROOT) $(LIBSROOT) 
Islands.o: Islands.hpp Islands.cpp
	$(CC) $(CFLAGS) -c Islands.cpp
clean:
	rm *.o 