===
CFD
===

.. contents::
   :local:

.. highlight:: latex

What is Fluid Mechanics?
========================

A science and branch of physics concerning:

* **Fluids** - liquids or gases
* **Mechanics** - application of laws of force and motion

There are two branches of fluid mechanics:

* **Fluid statics** - fluids at rest
* **Fluid dynamics** - fluids in motion

What is Thermodynamics?
=======================

A science and branch of physics concerning:

* **Heat**
* **Work**

And their relation to variables such as:

* **Internal energy**
* **Enthalpy**
* **Entropy**

What is CFD?
============

A science and branch of fluid mechanics and thermodynamics to predict fluid flow using:

* Conservation laws
* Numerical methods
* Algorithms
* Digital computers

What are the reasons for using CFD?
===================================

Allows virtual experiments that in the real world would be:

* Difficult
* Dangerous
* Expensive
* Impossible

What are the steps in a CFD process?
====================================

* **Physical model** - describe physical model
* **Mathematical model** - describe equations that correspond with physical model
* **Numerical methods** - describe numerical methods that correspond with physical model and **implement in code**
* **Geometry/Grid** - describe grid that corresponds with physical model and **implement in code**
* **Numerical solver** - describe solver that corresponds with physical model and **implement in code** and **compute solution**
* **Verification** - establish solution validity
* **Validation** - compare with experimental data 

Why do vector calculus and linear algebra play important roles in CFD?
==========================================================================

Vector calculus is important in CFD because it allows **description** e.g.

* Vector field - velocity of a flow
* Divergence of a flow field - expansion or compression of a flow
* Curl of a flow field - rotation of a flow

Linear algebra is important in CFD because it allows **solution** of the description of the flow e.g.

* Eigenvalues - charactersitic speeds of a system of hyperbolic equations
* Eigenvectors - characteristic directions of a system of hyperbolic equations
* TDMA - solution of tri-diagonal systems (e.g. Navier-Stokes momentum equation)


What can go wrong in CFD?
=========================

* **Over-simplicifcation** - Simplify domain to look only at part of system not the whole system
* **Numerical error** - Numerics introduces dissipation or dispersion
* **Geometry/grid** - Errors due to coordinate transformation


