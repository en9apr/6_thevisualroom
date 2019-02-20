=====================
 The Riemann Problem
=====================

The **Riemann Problem** is a family of problems, (which includes the 1D Euler Equations). 

There are two components for the Riemann problem:

1) A conservation law (e.g. the 1D Euler Equations)
2) Piecewise constant initial data with a single jump discontinuity

e.g. for Euler Equations:

For :math:`x \lt x_0`  :math:`\quad \mathbf{U}(x, t_0) = \mathbf{U}_L`

For :math:`x \gt x_0`  :math:`\quad \mathbf{U}(x, t_0) = \mathbf{U}_R`

:math:`\mathbf{U}_L` and :math:`\mathbf{U}_R` are constant vectors

Take :math:`x_0 = 0` and :math:`t_0=0`

Usefulness of the Riemann Problem
=================================

It has an exact analytical solution for the Euler Equations (and any scalar conservation laws, or any linear system of equations)

Features of the Riemann Problem
===============================

1) Solution is **self-similar** - solution stretches in space and time, but doesn't change shape.

* :math:`\mathbf{U}(x, t_1)` and :math:`\mathbf{U}(x, t_2)` are "similar"

* In other words, it really depends on a single variable :math:`x / t`

* Or the solution is constant along a characteristic line :math:`x = ct` (c = constant) on x-t plane

* Generally, self-similarity means only **one** independent variable like: :math:`x / t` or :math:`t / \sqrt{x}`, similar to the Blasius solution from boundary layer theory.

2) There is an analytical solution :math:`\Rightarrow` Riemann problem useful to test numerical schemes

3) Riemann problem appears as part of numerical formulation of several CFD methods, e.g. "wave capturing" methods :math:`\Rightarrow` "Riemann solvers" are at the heart of the method - may need to solve these multiple times - major cost of the solution

