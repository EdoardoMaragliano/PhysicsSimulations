CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall 
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`



frattale.x: frattale.cpp ReticoloF.o
	$(CC) $(CFLAGS) -o frattale.x frattale.cpp ReticoloF.o $(CFLAGSROOT) $(LIBSROOT)
ReticoloF.o: ReticoloF.hpp ReticoloF.cpp
	$(CC) $(CFLAGS) -c ReticoloF.cpp $(CFLAGSROOT) $(LIBSROOT) 
Islands.o: Islands.hpp Islands.cpp
	$(CC) $(CFLAGS) -c Islands.cpp
clean:
	rm *.o 