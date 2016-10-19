===========================
 Strong and Weak Solutions
===========================

.. contents::
   :local:

Conservative and Non-conservative Form of the Burgers' Equation
===============================================================

Non-Conservative Form
---------------------

.. math:: {\partial u \over \partial t} + u {\partial u \over \partial x} = 0


Conservative Form
-----------------

.. math:: {\partial u \over \partial t} + {\partial F \over \partial x} = 0

where the Flux is :math:`F` is (for Burgers' Equation): 

.. math:: F = {u^2 \over 2}

May also write:

.. math:: {\partial u \over \partial t} + A {\partial u \over \partial x} = 0

where the Jacobian :math:`A` is (for Burgers' Equation):

.. math:: A = A(u) = u

For more than 1 dimension, A is a matrix

.. math:: A = {\partial F_i \over \partial u_j}

Hyperbolic
==========

The Equation is **hyperbolic** :math:`\Rightarrow` all of the eigenvalues of A are real.

Strong Solution
===============

u is continuous, but bounded discontinuties in derivatives may occur

Weak Solution
=============

u is continuous almost everywhere, i.e. it has a jump (e.g. a shock wave)

