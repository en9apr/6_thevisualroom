# $MYHPSC/codes/mpi/Makefile

# Usage:  to try out test1.f90, for example:
#     make test FNAME=test1
# or
#     make test FNAME=test1 NUM_PROCS=8

NUM_PROCS ?= 4   # default if not specified on command line or env variable

FC = mpif90
FFLAGS = 
LFLAGS = 

.PHONY: clean, test

FNAME ?= test

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

$(FNAME).exe: $(FNAME).f90
	$(FC) $(FFLAGS) $(LFLAGS) $(FNAME).f90 -o $(FNAME).exe

test: $(FNAME).exe
	mpiexec -n $(NUM_PROCS) ./$(FNAME).exe

clean:
	rm -f *.o *.exe
