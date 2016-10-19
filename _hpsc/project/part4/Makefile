# $MYHPSC/project/part4/Makefile

OBJECTS = random_util.o problem_description.o mc_walk.o laplace_mc.o
OBJECTS2 = random_util.o problem_description.o mc_walk.o test1.o
OBJECTS2b = random_util.o problem_description.o mc_walk.o test1b.o
OBJECTS3 = random_util.o problem_description.o mc_walk.o test2.o
NUMPROCS ?= 4
FLAGS ?= 

.PHONY: mpi laplace test1 test1b test2

mpi: testquad.exe
	mpiexec -n $(NUMPROCS) ./testquad.exe

test1: test1.exe
	./test1.exe

test1b: test1b.exe
	./test1b.exe

test2: test2.exe
	./test2.exe

testquad.exe: $(OBJECTS)
	mpif90 $(FLAGS) $(OBJECTS) -o testquad.exe

test1.exe: $(OBJECTS2)
	gfortran $(FFLAGS) $(OBJECTS2) -o test1.exe

test1b.exe: $(OBJECTS2b)
	gfortran $(FFLAGS) $(OBJECTS2b) -o test1b.exe

test2.exe: $(OBJECTS3)
	gfortran $(FFLAGS) $(OBJECTS3) -o test2.exe

%.o : %.f90
	mpif90 -c $(FLAGS) $< 

laplace: laplace.exe
	./laplace.exe

laplace.exe: $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o laplace.exe

#%.o : %.f90
#	gfortran $(FFLAGS) -c  $<

clean:
	rm -f *.o *.exe *.mod
