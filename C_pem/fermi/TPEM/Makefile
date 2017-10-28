LDFLAGS =  -L. -lm -ltpem
EXENAME = skyrmionlattice
CPPFLAGS = -O2 -Wall -DNO_LAPACK -DNO_GSL
CC = icc
fort = ifort

all: $(EXENAME) 

$(EXENAME): func.o libtpem.a skyrmionlattice.o
	$(fort) -w skyrmionlattice.o func.o libtpem.a -o ../$(EXENAME) 


clean:
	rm -f core* $(EXENAME) *.o

libtpem.a: tpem.o tpem_sparse.o tpem_subspace.o
	rm -f libtpem.a
	ar rv libtpem.a tpem.o tpem_sparse.o tpem_subspace.o
	ranlib libtpem.a
	
tpem.o: tpem.c tpem.h
	$(CC) $(CPPFLAGS) -c tpem.c

tpem_sparse.o: tpem_sparse.c tpem.h
	$(CC) $(CPPFLAGS) -c tpem_sparse.c

tpem_subspace.o: tpem_subspace.c tpem.h
	$(CC) $(CPPFLAGS) -c tpem_subspace.c
func.o: func.c tpem.h 
	$(CC) $(CPPFLAGS) -c func.c
skyrmionlattice.o: skyrmionlattice.f90
	$(fort) -w -c skyrmionlattice.f90
