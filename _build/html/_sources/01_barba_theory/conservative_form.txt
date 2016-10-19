=================
Conservative Form
=================

.. contents::
   :local:

Conversion from Non-Conservative Form to Conservative Form
==========================================================

* Scalar conservation laws have the following **Conservative Form** in 1D (where :math:`F` = Flux):

.. math:: {{\partial u} \over {\partial t}} + {\partial \over {\partial x}} \left[ F(u(x,t))  \right] = 0 

* Scalar conservation laws have the following **Non-Conservative Form** in 1D (where :math:`A` = Jacobian):

.. math:: {{\partial u} \over {\partial t}} + A(u(x,t)) {{\partial u} \over {\partial x}}= 0 

* To find the Conservative Form from the Non-Conservative Form, we find the Flux, by equating the appropriate terms from the Conservative and Non-Conservative Forms:

.. math:: {\partial \over {\partial x}} (F(u)) = A(u) {{\partial u} \over {\partial x}}

* LHS given by the definition of a function of a function:

.. math:: {\partial \over {\partial x}} (F(u(x,t))) = 
   {{\partial F} \over {\partial u}} {{\partial u} \over {\partial x}} = 
   A(u) {{\partial u} \over {\partial x}}

* By comparing the above, we observe that:

.. math:: {{\partial F} \over {\partial u}} = A(u)

* By integration:

.. math:: F(u) = \int A(u) du

Conversion of the Burgers Equation from Non-Conservative to Conservative Form
=============================================================================

* Non-Conservative Form:

.. math:: {\partial u \over \partial t} + u {\partial u \over \partial x} = 0

* Identification of A(u):

.. math:: A(u) = u 

* Integration of A(u) to determine the flux F(u) 

.. math:: F(u) = \int A(u) du = \int u du = {{u^2} \over 2}

* Conservative Form:

.. math:: {\partial u \over \partial t} + {\partial \over {\partial x}}{ \left( {u^2 \over 2} \right)} = 0


**Conversion from Non-Conservative Form to Conservative Form could result in multiple Conservative Forms e.g.**

.. math:: u {\partial \over {\partial t}} + u^2 {\partial \over {\partial x}} = 0 

Converts to

.. math:: {\partial \over {\partial t}} \left( u^2 \over 2 \right) + 
   {\partial \over {\partial x}}{ \left( {u^3 \over 3} \right)} = 0 

**Is this now an alternate solution?**

What is the Difference Between the Conservative and Non-Conservative Forms?
===========================================================================

Using BD for the Flux term for both Conservative and Non-Conservative Forms (we had to do this as it doesn't make sense to integrate the Non-Conservative Form):

Conservative Form
-----------------

.. math:: \sum {\partial \over {\partial x}} \left( {{u^2} \over 2} \right) = 
   {{u_1^2-u_0^2} \over {2 \Delta x}} + 
   {{u_2^2-u_1^2} \over {2 \Delta x}} +
   {{u_3^2-u_2^2} \over {2 \Delta x}} = 
   {{u_3^2-u_0^2} \over {2 \Delta x}}

i.e. Conservative Form gives :math:`\text{flux}_{out} - \text{flux}_{in}` which is how it should be

In other words, if:

.. math:: {\partial u \over \partial t} + {\partial \over \partial x} (F(u)) = 0

Then the time derivative of the integral of the **conserved quantity** equals the difference between the outlet and inlet fluxes.

.. math:: {d \over dt} \int_{out}^{in} u(x,t) dx = F(u(out,t))-F(u(in,t))

Non-Conservative Form
---------------------

.. math:: \sum u {{\partial u} \over {\partial x}}= 
   u_1 {(u_1-u_0) \over {\Delta x}} + 
   u_2 {(u_2-u_1) \over {\Delta x}} +
   u_3 {(u_3-u_2) \over {\Delta x}} = 
   {{u_1 (u_1-u_0) + u_2 (u_2-u_1) + u_3 (u_3-u_2)} \over {\Delta x}}

i.e. Non-Conservative Form gives a growing function that is not conserving the flux! These are internal numerical sources that do not cancel out like in the conservative form.
