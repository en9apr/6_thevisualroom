Repository for the University of Washington's course on High Performance Computing.

The branch master contains the following files:

* **homework 1:**
     * **test1.py** (python script to test for environment variables and python version)
     * **test2.sh** (shell script to test for environment variables, gfortran and ipython version)
     * **test3.f90** (a dummy fortran code)

* **homework2:**
      * **hw2a.py** (python script for quadratic interpolation)
      * **hw2b.py** (python script for quadratic, cubic and polynomial interpolation)

* **homework 3:**
      * **newton.py** (ipython module to compute square root of a number using Newton's method)
      * **intersections.py** (ipython script to plot a graph where two functions intersect - uses newton.py)
      * **test1.f90, newton.f90, intersections.f90, functions.f90, Makefile** (fortran code for the intersection problem using max iterations for convergence)
      * **problem7/newton.f90, problem7/functions.f90, problem7/test_quartic.f90, problem7/Makefile, problem7/intersections.py** (intersection problem using residual error for convergence)

* **homework 4:**
      * **quadrature.f90** (module for integrating a function)
      * **test1.f90** (cubic function - uses quadrature module)
      * **test2.f90** (sinusoidal function - uses quadrature module)
      * **quadrature_omp.f90, test2_omp.f90** (same as above but using OpenMP)

* **homework 5:**
      * **Makefile** (makefile for the 1D solution)
      * **functions.f90** (sinusoial function module)
      * **quadrature.f90, test.f90** (trapezoid rule module + test - serial)
      * **quadrature2.f90, test2.f90** (simpson's rule module + test - serial)
      * **quadrature3.f90, test3.f90** (trapezoid rule module + test - OpenMP)
      * **quad2d/functions.f90, quad2d/quadrature.f90, quad2d/test.f90, quad2d/Makefile** (integrate a sinusoid with trapezoid rule in 2D with Open MP and load balancing)

* **homework 6:**
      * **Makefile** (makefile for this directory)
      * **copyvalue.f90** (passes a value from Process 0 to Processes 1-3 using MPI)
      * **copyvalue2.f90** (passes a value from Process 3 to Processes 0-2 using MPI)
      * **part2/functions.f90, part2/quadrature.f90, part2/Makefile** (the sinusoidal function, quadrature module and makefile)
      * **part2/test.f90** (compute the same integral over all Processes using MPI)
      * **part2/test2.f90** (integral using Master-Worker paradigm)

* **project/part 1:**
      * **functions.f90, quadrature.f90, test2.f90, Makefile** (compute a 1D integral using MPI by dividing up intervals between an indeterminate number of processes)
* **project/part 2:**
      * **functions.f90, quadrature_mc.f90, random_util.f90, test_quad_mc.f90, Makefile** (compute integral using Monte Carlo method for a 20 dimensional integral - serial)
      * **plot_mc_quad_error.py** (plot the result)
* **project/part 3:**
      * **laplace_mc.py** (python version)
      * **problem_description.f90, laplace_mc.f90, mc_walk.f90, random_util.f90, Makefile** (compute the solution to Laplace's Equation using the Monte Carlo Method - serial)
      * **test1.f90, test2.f90** (trials I made along the way)
      * **plot_mc_quad_error.py** (plot the result)
* **project/part 4:**
      * **problem_description.f90, laplace_mc.f90, mc_walk.f90, random_util.f90, Makefile** (compute the solution to Laplace's Equation using the Monte Carlo Method - using MPI)
      * **test1.f90, test1b.f90, test2.f90** (trials I made along the way)
      * **plot_mc_quad_error.py** (plot the result)