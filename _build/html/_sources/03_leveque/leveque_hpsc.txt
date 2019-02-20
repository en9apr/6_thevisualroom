=========
HPSC MOOC
=========

These projects are based on the course "High Performance Scientific Computing" by `Professor Randall LeVeque <https://amath.washington.edu/people/randy-leveque/>`_.

This page contains the following codes:

* `Workspace 1: Testing Scripts and Dummy Code`_
     * **test1.py** (python script to test for environment variables and python version)
     * **test2.sh** (shell script to test for environment variables, gfortran and ipython version)
     * **test3.f90** (a dummy fortran code)

* `Workspace 2: Quadratic, Cubic and Polynomial Interpolation`_
      * **hw2a.py** (python script for quadratic interpolation)
      * **hw2b.py** (python script for quadratic, cubic and polynomial interpolation)

* `Workspace 3: Newton's Method and Convergence Criteria`_
      * **newton.py** (ipython module to compute square root of a number using Newton's method)
      * **intersections.py** (ipython script to plot a graph where two functions intersect - uses newton.py)
      * **test1.f90, newton.f90, intersections.f90, functions.f90, Makefile** (fortran code for the intersection problem using max iterations for convergence)
      * **problem7/newton.f90, problem7/functions.f90, problem7/test_quartic.f90, problem7/Makefile, problem7/intersections.py** (intersection problem using residual error for convergence)

* `Workspace 4: OpenMP`_ 
      * **quadrature.f90** (module for integrating a function)
      * **test1.f90** (cubic function - uses quadrature module)
      * **test2.f90** (sinusoidal function - uses quadrature module)
      * **quadrature_omp.f90, test2_omp.f90** (same as above but using OpenMP)

* `Workspace 5: Load Balancing`_
      * **Makefile** (makefile for the 1D solution)
      * **functions.f90** (sinusoidal function module)
      * **quadrature.f90, test.f90** (trapezoid rule module + test - serial)
      * **quadrature2.f90, test2.f90** (simpson's rule module + test - serial)
      * **quadrature3.f90, test3.f90** (trapezoid rule module + test - OpenMP)
      * **quad2d/functions.f90, quad2d/quadrature.f90, quad2d/test.f90, quad2d/Makefile** (integrate a sinusoid with trapezoid rule in 2D with Open MP and load balancing)

* `Workspace 6: MPI`_
      * **Makefile** (makefile for this directory)
      * **copyvalue.f90** (passes a value from Process 0 to Processes 1-3 using MPI)
      * **copyvalue2.f90** (passes a value from Process 3 to Processes 0-2 using MPI)
      * **part2/functions.f90, part2/quadrature.f90, part2/Makefile** (the sinusoidal function, quadrature module and makefile)
      * **part2/test.f90** (compute the same integral over all Processes using MPI)
      * **part2/test2.f90** (integral using Master-Worker paradigm)

* `Project 1: 1D Integral and MPI`_
      * **functions.f90, quadrature.f90, test2.f90, Makefile** (compute a 1D integral using MPI by dividing up intervals between an indeterminate number of processes)

* `Project 2: 20D Integral and the Monte Carlo Method`_
      * **functions.f90, quadrature_mc.f90, random_util.f90, test_quad_mc.f90, Makefile** (compute integral using Monte Carlo method for a 20 dimensional integral - serial)
      * **plot_mc_quad_error.py** (plot the result)

* `Project 3: Laplace's Equation and the Monte Carlo Method`_
      * **laplace_mc.py** (python version)
      * **problem_description.f90, laplace_mc.f90, mc_walk.f90, random_util.f90, Makefile** (compute the solution to Laplace's Equation using the Monte Carlo Method - serial)
      * **test1.f90, test2.f90** (trials I made along the way)
      * **plot_mc_quad_error.py** (plot the result)
* `Project 4: Laplace's Equation, MPI and the Monte Carlo Method`_
      * **problem_description.f90, laplace_mc.f90, mc_walk.f90, random_util.f90, Makefile** (compute the solution to Laplace's Equation using the Monte Carlo Method - using MPI)
      * **test1.f90, test1b.f90, test2.f90** (trials I made along the way)
      * **plot_mc_quad_error.py** (plot the result)

Workspace 1: Testing Scripts and Dummy Code
===========================================
* **test1.py** (python script to test for environment variables and python version)

.. literalinclude:: ../_hpsc/homework1/test1.py

* **test2.sh** (shell script to test for environment variables, gfortran and ipython version)

.. literalinclude:: ../_hpsc/homework1/test2.sh

* **test3.f90** (a dummy fortran code)

.. literalinclude:: ../_hpsc/homework1/test3.f90
   :language: fortran

Workspace 2: Quadratic, Cubic and Polynomial Interpolation
==========================================================
* **hw2a.py** (python script for quadratic interpolation)

.. literalinclude:: ../_hpsc/homework2/hw2a.py

* **hw2b.py** (python script for quadratic, cubic and polynomial interpolation)

.. literalinclude:: ../_hpsc/homework2/hw2b.py

Workspace 3: Newton's Method and Convergence Criteria
=====================================================

* **newton.py** (ipython module to compute square root of a number using Newton's method)

.. literalinclude:: ../_hpsc/homework3/newton.py

* **intersections.py** (ipython script to plot a graph where two functions intersect - uses newton.py)

.. literalinclude:: ../_hpsc/homework3/intersections.py

* **test1.f90, newton.f90, intersections.f90, functions.f90, Makefile** (fortran code for the intersection problem using max iterations for convergence)

.. literalinclude:: ../_hpsc/homework3/test1.f90
   :language: fortran
   
.. literalinclude:: ../_hpsc/homework3/newton.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework3/intersections.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework3/functions.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework3/Makefile

* **problem7/newton.f90, problem7/functions.f90, problem7/test_quartic.f90, problem7/Makefile, problem7/intersections.py** (intersection problem using residual error for convergence)

.. literalinclude:: ../_hpsc/homework3/problem7/newton.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework3/problem7/functions.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework3/problem7/test_quartic.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework3/problem7/Makefile
 
.. literalinclude:: ../_hpsc/homework3/problem7/intersections.py

Workspace 4: OpenMP
===================

* **quadrature.f90** (module for integrating a function)

.. literalinclude:: ../_hpsc/homework4/quadrature.f90
   :language: fortran

* **test1.f90** (cubic function - uses quadrature module)

.. literalinclude:: ../_hpsc/homework4/test1.f90
   :language: fortran

* **test2.f90** (sinusoidal function - uses quadrature module)
 
.. literalinclude:: ../_hpsc/homework4/test2.f90
   :language: fortran

* **quadrature_omp.f90, test2_omp.f90** (same as above but using OpenMP)

.. literalinclude:: ../_hpsc/homework4/quadrature_omp.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework4/test2_omp.f90
   :language: fortran

Workspace 5: Load Balancing
===========================

* **Makefile** (makefile for the 1D solution)

.. literalinclude:: ../_hpsc/homework5/Makefile

* **functions.f90** (sinusoial function module)

functions.f90

.. literalinclude:: ../_hpsc/homework5/functions.f90
   :language: fortran

* **quadrature.f90, test.f90** (trapezoid rule module + test - serial)

.. literalinclude:: ../_hpsc/homework5/quadrature.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework5/test.f90
   :language: fortran

* **quadrature2.f90, test2.f90** (simpson's rule module + test - serial)

.. literalinclude:: ../_hpsc/homework5/quadrature2.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework5/test2.f90
   :language: fortran

* **quadrature3.f90, test3.f90** (trapezoid rule module + test - OpenMP)
 
.. literalinclude:: ../_hpsc/homework5/quadrature3.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework5/test3.f90
   :language: fortran

* **quad2d/functions.f90, quad2d/quadrature.f90, quad2d/test.f90, quad2d/Makefile** (integrate a sinusoid with trapezoid rule in 2D with Open MP and load balancing)

.. literalinclude:: ../_hpsc/homework5/quad2d/functions.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework5/quad2d/quadrature.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework5/quad2d/test.f90
   :language: fortran

Workspace 6: MPI
================

* **Makefile** (makefile for this directory)

.. literalinclude:: ../_hpsc/homework6/Makefile

* **copyvalue.f90** (passes a value from Process 0 to Processes 1-3 using MPI)

.. literalinclude:: ../_hpsc/homework6/copyvalue.f90
   :language: fortran

* **copyvalue2.f90** (passes a value from Process 3 to Processes 0-2 using MPI)

.. literalinclude:: ../_hpsc/homework6/copyvalue2.f90
   :language: fortran

* **part2/functions.f90, part2/quadrature.f90, part2/Makefile** (the sinusoidal function, quadrature module and makefile)

.. literalinclude:: ../_hpsc/homework6/part2/functions.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework6/part2/quadrature.f90
   :language: fortran

.. literalinclude:: ../_hpsc/homework6/part2/Makefile

* **part2/test.f90** (compute the same integral over all Processes using MPI)
 
.. literalinclude:: ../_hpsc/homework6/part2/test.f90
   :language: fortran

* **part2/test2.f90** (integral using Master-Worker paradigm)Makefile for this directory

.. literalinclude:: ../_hpsc/homework6/part2/test2.f90
   :language: fortran

Project 1: 1D Integral and MPI
==============================

* **functions.f90, quadrature.f90, test2.f90, Makefile** (compute a 1D integral using MPI by dividing up intervals between an indeterminate number of processes)

.. literalinclude:: ../_hpsc/project/part1/Makefile

.. literalinclude:: ../_hpsc/project/part1/functions.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part1/quadrature.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part1/test2.f90
   :language: fortran

Project 2: 20D Integral and the Monte Carlo Method
==================================================

* **functions.f90, quadrature_mc.f90, random_util.f90, test_quad_mc.f90, Makefile** (compute integral using Monte Carlo method for a 20 dimensional integral - serial)
* **plot_mc_quad_error.py** (plot the result)

.. literalinclude:: ../_hpsc/project/part2/Makefile

.. literalinclude:: ../_hpsc/project/part2/functions.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part2/quadrature_mc.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part2/random_util.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part2/test_quad_mc.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part2/plot_mc_quad_error.py

Project 3: Laplace's Equation and the Monte Carlo Method
========================================================

* **laplace_mc.py** (python version)
* **problem_description.f90, laplace_mc.f90, mc_walk.f90, random_util.f90, Makefile** (compute the solution to Laplace's Equation using the Monte Carlo Method - serial)
* **test1.f90, test2.f90** (trials I made along the way)
* **plot_mc_quad_error.py** (plot the result)

.. literalinclude:: ../_hpsc/project/part3/laplace_mc.py

.. literalinclude:: ../_hpsc/project/part3/Makefile

.. literalinclude:: ../_hpsc/project/part3/problem_description.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part3/laplace_mc.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part3/mc_walk.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part3/random_util.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part3/test1.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part3/test2.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part3/plot_mc_quad_error.py

Project 4: Laplace's Equation, MPI and the Monte Carlo Method
=============================================================

* **problem_description.f90, laplace_mc.f90, mc_walk.f90, random_util.f90, Makefile** (compute the solution to Laplace's Equation using the Monte Carlo Method - using MPI)
* **test1.f90, test1b.f90, test2.f90** (trials I made along the way)
* **plot_mc_quad_error.py** (plot the result)

.. literalinclude:: ../_hpsc/project/part4/Makefile

.. literalinclude:: ../_hpsc/project/part4/problem_description.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part4/laplace_mc.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part4/mc_walk.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part4/random_util.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part4/test1.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part4/test1b.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part4/test2.f90
   :language: fortran

.. literalinclude:: ../_hpsc/project/part4/plot_mc_quad_error.py
