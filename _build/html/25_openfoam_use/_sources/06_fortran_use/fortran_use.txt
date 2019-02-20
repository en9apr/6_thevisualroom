===============
 Using Fortran
===============

.. contents::
   :local:

.. highlight:: fortran

Developing
==========

1) Understand the Problem

   - Question being asked
   - Mathematical description

2) Formulate the Problem

   - Inputs
   - Outputs
   - Variables and constants
   - Flow chart

3) Design Algorithm to Solve the Problem

   - Numerical scheme
   - Pseudo-code

4) Implement Algorithm in Code

**Steps 3) and 4) may be incremental - where each sub-step can be tested**

Useful to use comments at each step so you are clear about the question being asked and what you are going to do (like a quick step 1) and 2))

Always use meaningful names for programs, modules, functions, subroutines and variables.

Naming Convections
==================

1) Use lowercase for all Fortran constructs (do, subroutine, module, ...).
2) Follow short mathematical notation for mathematical variables (sigma, gamma, rho, epsilon, ...).
3) For other names use all lowercase: try to keep names to one or two syllables; if more are required, use underscores to clarify (spline_interpolate, stop_error).

This convection is quite good for Fortran 90:

`Fortran Naming Convection <http://www.fortran90.org/src/best-practices.html>`_

Program Stucture
----------------

* The trick with names is to avoid ambigous interpretation, e.g. avoid:

 - ``mod`` could be module or modulus
 - ``int`` could mean integrate or interpolate
 - ``sub`` could be subroutine or subset
 - ``pro`` could be program or product

* Fortran usually prefers whole words, but these may become too verbose at some point

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           program area_under_curve

     - Program names usually denote the overall task - "what is being tested?" usually answers this one. Prefixes like ``test_`` can sometimes be useful for programs
   * - ::

           module integrate
           module problem_description

     - Module names should describe the overall function of the procedures within that module
   * - ::

           function trapezoid_rule(x, y ,z)
           function simpsons_rule(x, y, z)
           
     - Function names are usually verb-noun combinations. Prefixes like ``print_``  ``write_``  ``get_`` and ``set_`` are sometimes useful for functions.
   * - ::

           subroutine initial_conditions(u, v, x, y, t)
           subroutine analytical_solution(u, v, x, y, t)
           subroutine boundary_conditions(u, v, x, y, t)
           
     - Similar convention for subroutines as for functions.


Variables and Constants
-----------------------

As these are used in formulae, it's important to keep them as short as possible. Note that Fortran is **case insensitive**

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

            x, y, z 

     - Cartesian coordinates
   * - ::
     
            i, j, k 

     - Usually denote spatial coordinates
   * - ::
     
            n, m
 
     - Usually denote time or interation number
   * - ::

            process_no, step_no

     - ``_no`` can denote number (as is usual in english)
   * - ::

            i_max, j_max, k_max, n_max, m_max, process_no_max, step_no_max
     - Maximum value or range denoted by ``_max``. Max denoting maximum as is usual
   * - ::

            PI, RHO, NU

     - Constants (parameters) can be in uppercase


Indentation
===========

Use 4 spaces for indentation

Templates
=========

Module Template
---------------

* ``implicit none`` works for whole module
* ``private`` implies that everything in this module is private by default
* Only make public what you want

::

    module integrate

        use constants, only: RHO, PI
        use utilities, only: stop_error
        implicit none
        private
        public integrate, normalize, &
               euler, runge_kutta

    contains

        subroutine get_values(a, b, c)
    
        end subroutine
    
    end module

Program Template
----------------

* Note the use of **explicit** imports, i.e. avoid ``use integrate``. Instead say what you are using ``use integrate, only: midpoint``

::

    program uranium
        use mesh, only: create_mesh
        use utilities, only: stop_error, RHO
        implicit none

        integer, parameter :: Z = 92
        real(kind=8), parameter :: R_MIN = 8.0d-9, R_MAX = 50.0d-10, A = 1.0d+7
    
        print *, "I am running"
    end program

Comments
========

Comments have two functions - explain the interface and explain the implementation

Interface
---------

Explain what is does - the **interface**. Put it before any code is written

Functions/Subroutines
~~~~~~~~~~~~~~~~~~~~~

::

    ! Estimate the integral of f(x) from a to b using the
    ! Trapezoid Rule with n points.

    ! Input:
    !   f:  the function to integrate
    !   a:  left endpoint
    !   b:  right endpoint
    !   n:  number of points to use
    ! Output:
    !   the estimate of the integral

    function trapezoid


Modules
~~~~~~~

::

    ! Performs integration using quadrature integration 
    ! and creates a table of the error between this and
    ! the known solution.

    module quadrature_omp  


Programs
~~~~~~~~

::

    ! Prints a table for the effect of the number
    ! of integration points on the accuracy of the integration

    ! Example use: 
    ! $ gfortran -fopenmp quadrature_omp.f90 test2_omp.f90
    ! $ ./a.out

    program test_openmp


Implementation
--------------

Explain how it does it - the **implementation**

Above a block of code to denote how it does it (never usually need these inline - that's too  much commenting)

::
 
    ! Print the number of function evaluations by each thread:
        do i=0,nthreads-1
            print 101,  i, fevals(i)
    101     format("fevals by thread ",i2,": ",i13)
            enddo

Makefiles
=========

The general form of a make command is as follows. Use a tab character not spaces:

::

   target: dependencies
       <TAB> commands(s) to make target

Example where there is no output file. **It's important to know the dependencies for the module files, as they will compile in that order**, i.e. ``MODULES = functions.mod newton.mod`` implies that the newton module depends on the functions module.

::

    # $MYHPSC/homework3/Makefile
    # Example usage:
    #	$ make test1
    #	$ make clean

    # Dependencies for test1.exe - The names of object files (\ is a continuation character):
    OBJECTS = functions.o \
              newton.o \
              test1.o

    # Dependencies for test1.exe - The names of modules (needed if modules are used):
    MODULES = functions.mod \
              newton.mod

    # FLAGS
    # -c flag means compile to one file (very common if program is in many files)
    # -o FILENAME means rename output from a.out to FILENAME.exe
    # -g generates extra debugging information usable by GDB
    # -03 level 3 optimisation for compling code

    FFLAGS = -g

    # Phony targets don't create files (e.g deletes files or prints to screen)
    .PHONY: test1 clean

    # 1) Highest level: dependency for make test1 is test1.exe
    #    If older it runs ./test1.exe
    test1: test1.exe
    	./test1.exe

    # 2) Second highest level: dependencies for test1.exe are .mod and .o files
    #    If older it runs gfortran complier
    test1.exe: $(MODULES) $(OBJECTS)
    	gfortran $(FFLAGS) $(OBJECTS) -o test1.exe

    # 3) Third highest level: dependencies for .o and .mod files are .f90 files
    #    If older it runs gfortran complier ($< refers to dependency)
    %.o : %.f90
    	gfortran $(FFLAGS) -c  $< 

    %.mod: %.f90
    	gfortran $(FFLAGS) -c $<

    # Removes all files (rm) -f means force nonexistent files (never prompt)
    clean:
    	rm -f *.o *.exe *.mod


You can have multiple Makefiles in the same directory. To distinguish use:

::

    $ make test1 -f Makefile_2

To print to the screen from a Makefile use ``@echo``

::

    OBJECTS = functions.o newton.o test1.o
    MODULES = functions.mod newton.mod
    .PHONY: test

    test:
            @echo "Modules are: " $(MODULES)
            @echo "Objects are: " $(OBJECTS)

Compiling
=========

**If the Makefile does not compile, you may need to compile separately to see what the errors are in the compilation**

Simple compilation (single file)
--------------------------------

::

   $ gfortran program.f90

Multifile complitation (with modules)
-------------------------------------

* Compile first (separate compile needed because of module dependencies):

::

   $ gfortran -c program.f90 
   $ gfortran -c module_1.f90 
   $ gfortran -c module_2.f90

* Then link (with optional executable rename):

::

   $ gfortran program.o module_1.o module_2.o -o program.exe

* Easier to use Makefile for multi-file programs

Flags
=====

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           -c

     - Compile to an object file, object files are later linked into a complete program
   * - ::

           -o FILENAME

     - Specifies the name of the output file
   * - ::

           -fopenmp

     - Compile using OpenMP
   * - ::

           -g

     - Generates extra debugging information for GDB
   * - ::

           -g3

     - Generates even more debugging information
   * - ::

           -03

     - Optimised code - program is faster but takes longer to compile

Timing
======

Total CPU Time
--------------

This is for the total CPU time (called the "user" time). The real time and sys time aren't really important.

::

    $ time ./a.out
    <output from code>

    real    0m5.279s
    user    0m1.915s
    sys     0m0.006s

CPU Time for Part of Code
-------------------------

The ``cpu_time`` tells the CPU time used between two successive calls:

::

    real(kind=8) :: t1, t2, elapsed_cpu_time

    call cpu_time(t1)

    !code to be timed

    call cpu_time(t2)
    elapsed_cpu_time = t2 - t1

CPU Time is proportional to :math:`n^3` if using matrix multiplication, doubling n will multiply CPU time by 8 (although larger matrices may be affected by cache). Optimisation flags can be used for reducing CPU time.

Debugging
=========

Print Statements
----------------

Adding print statements to a program is a tried and true method of debugging, and the only method that many programmers use.

Print statements can be added almost anywhere in a Fortran code to print things out to the terminal window as it goes along.

You might want to put some special symbols in debugging statements to flag them as such, which makes it easier to see what output is your debug output, and also makes it easier to find them again later to remove from the code, e.g. you might use “+++” or “DEBUG”.

Compiling with  gfortran Flags
------------------------------

There are a number of flags you can use when compiling your code that will make it easier to debug.

Here’s a generic set of options you might try:

::

    $ gfortran -g -W -Wall -fbounds-check -pedantic-errors \
    -ffpe-trap=zero, invalid, overflow, underflow program.f90

Most of these options indicate that the program should give warnings or die if certain bad things happen.

Compiling with the -g flag indicates that information should be generated and saved during compilation that can be used to help debug the code using a debugger such as gdb or totalview. You generally have to compile with this option to use a debugger.

The gdb Debugger
----------------

gdb is an open source debugger, but often doesn't work well with Fortran. 

::

    $ cd $UWHPSC/codes/fortran 
    $ gfortran -g segfault1.f90 
    $ gdb a.out 

    (gdb) run

     Runs for a while and then prints
     "Program received signal EXC_BAD_ACCESS, 
     Could not access memory." Tells what line it died in.

     11 a(i) = 5
     (gdb) p i 
      $1 = 241 
     (gdb) q

This at least revels where the error happened and allows printing the value of i where it died.

Software
========

It is best to use high quality software as much as possible for several reasons:

1) It will take less time to figure out how to use the software then to write your own version (assuming it's well documented)

2) Good general software has been extensively tested on a wide variety of problems

3) Often general software is much more sophisticated than what you might write yourself, for example it might provide error estimates automatically, or it might be optimized to run fast.

Often we use:

* `LAPACK <http://www.netlib.org/lapack/>`_ (Linear Algebra Package)
* `BLAS <http://www.netlib.org/blas/>`_ (Basic Linear Algebra Subprograms)

There are others, such as:

* `Clawpack <http://www.clawpack.org/>`_. Clawpack stands for “Conservation Laws Package” and was initially developed for linear and nonlinear hyperbolic systems of conservation laws, with a focus on implementing high-resolution Godunov type methods using limiters in a general framework applicable to many applications. These finite volume methods require a “Riemann solver” to resolve the jump discontinuity at the interface between two grid cells into waves propagating into the neighboring cells. The formulation used in Clawpack allows easy extension to the solution of hyperbolic problems that are not in conservation form.

There are many version of Clawpack, such as:

* AMRClaw includes block-structured adaptive mesh refinement that allows one to use a non-uniform grid that changes in time and uses smaller grid cells in regions with fine structure or where high accuracy is required.
* GeoClaw Includes the AMR capabilities of AMRClaw and also has a number of special routines and algorithms for handling geophysical problems, including special well-balanced, positivity-preserving shallow water solvers.
* PyClaw includes the high-order WENO-RK algorithms of SharpClaw.
