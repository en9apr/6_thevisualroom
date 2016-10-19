
OBJECTS = functions.o quadrature.o test.o
OBJECTS2 = functions.o quadrature.o test2.o
FFLAGS = 
NUM_PROCS ?= 4   # default if not specified on command line or env variable

.PHONY: test clean test2

test: test.exe
	mpiexec -n $(NUM_PROCS) test.exe

test2: test2.exe
	mpiexec -n $(NUM_PROCS) test2.exe

test.exe: $(OBJECTS)
	mpif90 $(FFLAGS) $(OBJECTS) -o test.exe

test2.exe: $(OBJECTS2)
	mpif90 $(FFLAGS) $(OBJECTS2) -o test2.exe

%.o : %.f90
	mpif90 $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod

