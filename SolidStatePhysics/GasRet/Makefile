CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall 
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`



main.x: main.cpp Matrice.o
	$(CC) $(CFLAGS) -o main.x main.cpp Matrice.o $(CFLAGSROOT) $(LIBSROOT)
Matrice.o: Matrice.hpp Matrice.cpp
	$(CC) $(CFLAGS) -c Matrice.cpp $(CFLAGSROOT) $(LIBSROOT)
clean:
	rm *.o
clean-all:
	rm *.o *.x