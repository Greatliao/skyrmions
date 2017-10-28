fort = ifort
EXENAME = ../skyrmionlattice
CPPFLAGS = -r8 -w
LINK = -L/opt/intel/parallel/composer_xe_2013/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

all: $(EXENAME)

$(EXENAME):  skyrmionlattice.f90 nearest_site.o random_number.o
	$(fort) $(CPPFLAGS) skyrmionlattice.f90 nearest_site.o random_number.o -o $(EXENAME) $(LINK)


clean:
	rm -f core* $(EXENAME) *.o

nearest_site.o: nearest_site.f90
	$(fort) -c nearest_site.f90
random_number.o: random_number.f90
	$(fort) -c random_number.f90
