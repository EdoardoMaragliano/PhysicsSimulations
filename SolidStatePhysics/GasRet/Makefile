CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall 
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`



gas_reticolare.x: gas_reticolare.cpp Matrice.o
	$(CC) $(CFLAGS) -o gas_reticolare.x gas_reticolare.cpp Matrice.o $(CFLAGSROOT) $(LIBSROOT)
Matrice.o: Matrice.hpp Matrice.cpp
	$(CC) $(CFLAGS) -c Matrice.cpp $(CFLAGSROOT) $(LIBSROOT)

gret.x: gret.cpp
	$(CC) $(CFLAGS) -o gret.x gret.cpp $(CFLAGSROOT) $(LIBSROOT)

 

clean:
	rm *.o 