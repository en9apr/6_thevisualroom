==================================================================
Order of Accuracy, Midpoint Scheme and Model Equations
==================================================================

Order of Accuracy
=================

.. math:: \left .  {\partial u \over \partial x} \right \vert_i = \lim_{\Delta x \rightarrow 0} {u(x_i + \Delta x) - u(x_i) \over \Delta x}

* If :math:`\Delta x` is small :math:`\rightarrow` approximates :math:`\partial u / \partial x`
* Improve approximation if :math:`\Delta x` is reduced
* Errors are always introduced (truncation error)

Truncation Error
~~~~~~~~~~~~~~~~

* Definition: The power of :math:`\Delta x` with which the truncation error tends to zero is called the Order of Accuracy of the Finite Difference approximation.
* The Taylor Series Expansions: 

  - FD and BD are both first order or are :math:`O(\Delta x)` (Big-O Notation)
  - CD is second order or are :math:`O(\Delta x^2)` (Big-O Notation)

Difference formulas for first derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* FD and BD are called "one-sided" formulas, in other words they only involve points on one side of the point you are evaluating.
* FD:

.. math:: \left . {\partial u \over \partial x} \right \vert_i = {(u_{i+1} - u_i) \over \Delta x} -  {\Delta x \over 2} \left . {\partial^2 u \over \partial x^2} \right \vert_i - {\Delta x^2 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i - \cdots - h.o.t

* BD:

.. math:: \left . {\partial u \over \partial x} \right \vert_i = {(u_i - u_{i-1}) \over \Delta x} + {\Delta x \over 2} \left . {\partial^2 u \over \partial x^2} \right \vert_i - {\Delta x^2 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i - \cdots - h.o.t

* CD - Add FD & BD:

.. math:: \left . {\partial u \over \partial x} \right \vert_i = {(u_{i+1} - u_{i-1}) \over 2 \Delta x} - {\Delta x^2 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i - \cdots - h.o.t

Difference Formula for Second Derivative
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* FD:

.. math:: u_{i+1} = u_i + \Delta x \left . {\partial u \over \partial x} \right \vert_i + {\Delta x^2 \over 2} \left . {\partial^2 u \over \partial x^2} \right \vert_i + {\Delta x^3 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i + \cdots + h.o.t

* BD:

.. math:: u_{i-1} = u_i - \Delta x \left . {\partial u \over \partial x} \right \vert_i + {\Delta x^2 \over 2} \left . {\partial^2 u \over \partial x^2} \right \vert_i - {\Delta x^3 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i + \cdots + h.o.t

* CD - Add FD & BD:

.. math:: \left . {\partial^2 u \over \partial x^2} \right \vert_i = {(u_{i+1}-2u_i + u_{i-1}) \over  \Delta x^2} - O(\Delta x^2)

Meaning of the Accuracy
~~~~~~~~~~~~~~~~~~~~~~~

* Note: 1D domain (0,1) with 11 points (10 intervals) so :math:`\Delta x = 0.1`

- 1st order :math:`\sim O(\Delta x) \sim O(10 \%)` Error is order of 10%
- 2nd order :math:`\sim O(\Delta x ^2) \sim O(1 \%)` Error is order 1%
- For 1% in 1st order: :math:`\sim O(1 \%)` need 100 divisions :math:`\Delta x = 0.01` i.e. 101 mesh points
- Hence first order methods are more expensive than second order.

Midpoint Scheme
===============

* Now consider the CD approximation and look at a CD approximation at the point :math:`i+{1 \over 2}`

.. math:: \left . {\partial u \over \partial x} \right \vert_{i+{1 \over 2}} = {(u_{i+1} - u_{i}) \over \Delta x} - O(\Delta x^2)

* Now consider the CD approximation and look at a CD approximation at the point :math:`i-{1 \over 2}`

.. math:: \left . {\partial u \over \partial x} \right \vert_{i-{1 \over 2}} = {(u_{i} - u_{i-1}) \over \Delta x} - O(\Delta x^2)

* It looks like FD and BD gained an order of accuracy (they are really both just CD)

Model Equations
===============

* Recall the Navier-Stokes Equations:

.. math:: {\partial \vec V \over \partial t} + \vec V (\nabla \cdot \vec V) = -{\nabla p \over \rho} + \nu \nabla^2 \vec V

* Diffusive terms appear through 2nd order derivative terms. If inviscid :math:`\rightarrow 0`
* Convective fluxes appear as 1st order derivatives in space. If symmetry :math:`\rightarrow 0`
* Transient term appears as 1st order derivative. If steady :math:`\rightarrow 0`

* Various modelling assumptions (inviscid, symmetric, steady etc) :math:`\rightarrow` Model equations
* This results in a system of PDEs, where the highest space derivative is 2nd order and the highest time derivative is 1st order.


**How do we identify whether diffusion, convection and unsteadyness are compatible in terms of the physics?**

* Numerical discretizations: Must be adequate for the physical process :math:`\rightarrow` Convection, Diffusion, Unsteadyness

* Model Equations are simplified forms of the Navier Stokes Equations.
