CC=CC
COMPFLAGS= -c -Dlapack_complex_double="std::complex<double>"  -std=c++11 
DFLAGS =   -DTSQR -DDMD -DPrint
MPIFLAGS=
# MPIFLAGS=-lmpi++ -lmpi

# Possible directive flags: -DSingleProc -DOzgurTest -DOzgurN2Case  -DCovariancePOD -DPrint (-DDirectSVD or -DTSQR) -DCheckSVD -DCheckQR

EXECDIR = ../build
all: $(EXECDIR)/solver.x Makefile
	rm -rf *.o

 LIBLINK=-llapacke -llapack -lblas -lgfortran -lm $(MPIFLAGS)


$(EXECDIR)/solver.x: main.o Makefile
	$(CC) main.o -o $(EXECDIR)/solver_dmd.x $(LIBLINK)

main.o: main.cpp Makefile
	$(CC) $(COMPFLAGS) $(DFLAGS) main.cpp

clean:
	rm -rf *.o

veryclean:
	rm -rf *.o *.x
