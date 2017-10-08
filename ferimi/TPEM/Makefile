LDFLAGS =  -L. -lm -ltpem
EXENAME = tpem_test
CPPFLAGS = -O2 -Wall -DNO_LAPACK -DNO_GSL
CC = icc

all: $(EXENAME) 

$(EXENAME): tpem_test.o diag.o libtpem.a 
	$(CC) $(CPPFLAGS) tpem_test.o diag.o -o $(EXENAME) $(LDFLAGS) 

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
	
