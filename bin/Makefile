F90 = gfortran
CFLAGS = -O3 -fopenmp -Wpedantic
LDFLAGS = -O3 -fopenmp -Wpedantic
SRCDIR = ../src

OBJ = ni_types.o spherical_harmonics.o lebedev.o ni_fun.o ni_grid.o ni_module.o ni_gradients.o unit_test.o test_gradients.o test_suite.o

libni: $(OBJ)
	$(F90) $(LDFLAGS) $(OBJ) -o libni.exe

clean:
	rm -f *.o *.mod

%.o: $(SRCDIR)/%.f90
	$(F90) $(CFLAGS) -c $<
: