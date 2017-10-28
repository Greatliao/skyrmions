EXENAME = skyrmionlattice
fort = ifort

all: $(EXENAME) 

$(EXENAME): pem.o sparse_matrix.o nearest_site.o skyrmionlattice.o 
	$(fort) pem.o sparse_matrix.o nearest_site.o skyrmionlattice.o -o ../$(EXENAME)

clean:
	rm -f core* $(EXENAME) *.o

pem.o: pem.f90
	$(fort) -c pem.f90

sparse_matrix.o: sparse_matrix.f90
	$(fort) -c sparse_matrix.f90

nearest_site.o: nearest_site.f90
	$(fort) -c nearest_site.f90

skyrmionlattice.o: skyrmionlattice.f90
	$(fort) -w -c skyrmionlattice.f90

