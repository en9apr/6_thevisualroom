
OBJECTS = functions.o quadrature.o test.o
OBJECTS2 = functions.o quadrature.o test2.o
OBJECTS3 = functions.o quadrature.o test3.o
FFLAGS = 
NUM_PROCS ?= 4   # default if not specified on command line or env variable

.PHONY: test clean test2 test3

test: test.exe
	mpiexec -n $(NUM_PROCS) test.exe

test2: test2.exe
	mpiexec -n $(NUM_PROCS) test2.exe

test3: test3.exe
	mpiexec -n $(NUM_PROCS) test3.exe

test.exe: $(OBJECTS)
	mpif90 $(FFLAGS) $(OBJECTS) -o test.exe

test2.exe: $(OBJECTS2)
	mpif90 $(FFLAGS) $(OBJECTS2) -o test2.exe

test3.exe: $(OBJECTS3)
	mpif90 $(FFLAGS) $(OBJECTS3) -o test3.exe

%.o : %.f90
	mpif90 $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod

