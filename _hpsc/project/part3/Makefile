OBJECTS = random_util.o problem_description.o mc_walk.o test1.o
OBJECTS2 = random_util.o problem_description.o mc_walk.o test2.o
OBJECTS3 = random_util.o problem_description.o mc_walk.o laplace_mc.o
FFLAGS = 

.PHONY: test test2 laplace

test: test.exe
	./test.exe

test2: test2.exe
	./test2.exe

laplace: laplace.exe
	./laplace.exe

test.exe: $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test.exe

test2.exe: $(OBJECTS2)
	gfortran $(FFLAGS) $(OBJECTS2) -o test2.exe

laplace.exe: $(OBJECTS3)
	gfortran $(FFLAGS) $(OBJECTS3) -o laplace.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod

