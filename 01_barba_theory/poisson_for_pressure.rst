=============================
Poisson Equation for Pressure
=============================

For compressible flow, pressure and velocity can be coupled with the Equation of State. **But for incompressible flow, there is no obvious way to couple pressure and velocity.**

So, take the divergence of the momentum equation and use the continuity equation to get a Poisson equation for pressure.

Poisson's Equation in Continuous Domain
=======================================

Navier Stokes Equations in Vector Notation
------------------------------------------

The Continuity Equation for an **incompressible fluid**:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the divergence operator (div):

.. math:: \text{div}\ \mathbf{u} = 0

Or using the nabla operator :math:`\nabla`:

.. math:: \nabla \cdot \mathbf{u} = 0 
   :label: one

The Momentum Equation for an **incompressible fluid**:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the gradient (grad) and Laplacian operators (:math:`\nabla^2`):

.. math:: {{\partial \mathbf{u}} \over {\partial t}} + \mathbf{u} \cdot \text{grad } \mathbf{u} =
          {-{1 \over \rho} \text{grad } p} + {\nu \nabla^2 \mathbf{u}} 


Using the nabla operator :math:`\nabla` and Laplacian operators (:math:`\nabla^2`) :

.. math:: {{\partial \mathbf{u}} \over {\partial t}} + \mathbf{u} \cdot \nabla \mathbf{u} =
          {-{1 \over \rho} \nabla p} + {\nu \nabla^2 \mathbf{u}} 
   :label: two

* Mass conservation :eq:`one` for :math:`\rho = \text{constant}` is a kind of **kinematic constraint** to the momentum equation :eq:`two`

* PROBLEM: There is no obvious way to couple the velocity and pressure (for compressible fluids we would have an equation of state, which provides a relation between density and pressure)

Continuity Equation
~~~~~~~~~~~~~~~~~~~

.. math:: {{\partial u} \over {\partial x}} + {{\partial v} \over {\partial y}} = 0

Momentum Equations
~~~~~~~~~~~~~~~~~~

.. math:: {{\partial u} \over {\partial t}} + {u {{\partial u} \over {\partial x}}} + v {{\partial u} \over {\partial y}} = 
          -{1 \over \rho} {{\partial p} \over {\partial x}} +
          \nu \left( {{\partial^2 u} \over {\partial x^2}} + {{\partial^2 u} \over {\partial y^2}} \right)

.. math:: {{\partial v} \over {\partial t}} + {u {{\partial v} \over {\partial x}}} + v {{\partial v} \over {\partial y}} = 
          -{1 \over \rho} {{\partial p} \over {\partial y}} +
          \nu \left( {{\partial^2 v} \over {\partial x^2}} + {{\partial^2 v} \over {\partial y^2}} \right)

The Divergence of the Momentum Equation in Partial Notation
-----------------------------------------------------------

Now take the divergence of the momentum equations, if :math:`\mathbf{M}` is the "vector" of momentum equations, the divergence is:

.. math:: \nabla \cdot \mathbf{M} = {\partial \over {\partial x}} M_x + {\partial \over {\partial y}} M_y 


.. math:: {\partial \over {\partial x}} M_x =
          {\partial \over {\partial x}}  {{\partial u} \over {\partial t}} +
          {{\partial u} \over {\partial x}} {{\partial u} \over {\partial x}} + u {{\partial^2 u} \over {\partial x^2}} + 
          {{\partial v} \over {\partial x}} {{\partial u} \over {\partial y}} + v {{\partial^2 u} \over {\partial x \partial y}}  = 
          -{1 \over \rho} {{\partial^2 p} \over {\partial x^2}} +
          \nu \left( {{\partial^3 u} \over {\partial x^3}} + {{\partial^3 u} \over {\partial x \partial y^2}} \right)

.. math:: {\partial \over {\partial y}} M_y =
          {\partial \over {\partial y}}  {{\partial v} \over {\partial t}} +
          {{\partial u} \over {\partial y}} {{\partial v} \over {\partial x}} + u {{\partial^2 v} \over {\partial x \partial y}} + 
          {{\partial v} \over {\partial y}} {{\partial v} \over {\partial y}} + v {{\partial^2 y} \over {\partial y^2}}  = 
          -{1 \over \rho} {{\partial^2 p} \over {\partial y^2}} +
          \nu \left( {{\partial^3 v} \over {\partial x^2 \partial y}} + {{\partial^3 v} \over {\partial y^3}} \right)

Add the LHS
~~~~~~~~~~~

.. math:: {\partial \over {\partial x}} M_x + {\partial \over {\partial y}} M_y =
          {\partial \over {\partial t}} \left( {{\partial u} \over {\partial x}} + {{\partial v} \over {\partial y}} \right) +
          \left( {{\partial u} \over {\partial x}} \right)^2 + 
          {{\partial u} \over {\partial y}} {{\partial v} \over {\partial x}}+
          u {{\partial^2 u} \over {\partial x^2}} +
          u {{\partial^2 v} \over {\partial x \partial y}} +
          {{\partial v} \over {\partial x}} {{\partial u} \over {\partial y}} +
          \left( {{\partial v} \over {\partial y}} \right)^2 + 
          v {{\partial^2 u} \over {\partial x \partial y}} +
          v {{\partial^2 v} \over {\partial y^2}} = RHS

Re-arrange:

.. math:: {\partial \over {\partial x}} M_x + {\partial \over {\partial y}} M_y =
          {\partial \over {\partial t}} \left( {{\partial u} \over {\partial x}} + {{\partial v} \over {\partial y}} \right) +
          \left( {{\partial u} \over {\partial x}} \right)^2 + 
          2 {{\partial u} \over {\partial y}} {{\partial v} \over {\partial x}}+
          u {\partial \over {\partial x}} \left( {{\partial u} \over {\partial x}} +
          {{\partial v} \over {\partial y}} \right) +
          \left( {{\partial v} \over {\partial y}} \right)^2 + 
          v {\partial \over {\partial y}} \left( {{\partial u} \over {\partial x}} +
          {{\partial v} \over {\partial y}} \right) = RHS

Apply Continuity, so :math:`\nabla \cdot \mathbf{u} = 0`. Hence:

.. math:: {\partial \over {\partial x}} M_x + {\partial \over {\partial y}} M_y =
          \left( {{\partial u} \over {\partial x}} \right)^2 + 
          2 {{\partial u} \over {\partial y}} {{\partial v} \over {\partial x}}+
          \left( {{\partial v} \over {\partial y}} \right)^2 = RHS

Add the RHS
~~~~~~~~~~~

.. math:: -{1 \over \rho} \left( {{\partial^2 p} \over {\partial x^2}} + {{\partial^2 p} \over {\partial y^2}} \right)+
          \nu \left( {{\partial^3 u} \over {\partial x^3}} + {{\partial^3 u} \over {\partial x \partial y^2}} +
           {{\partial^3 v} \over {\partial x^2 \partial y}} + {{\partial^3 v} \over {\partial y^3}} \right) = LHS

Re-arrange:

.. math:: -{1 \over \rho} \left( {{\partial^2 p} \over {\partial x^2}} + {{\partial^2 p} \over {\partial y^2}} \right)+
          \nu \left( {{\partial^2} \over {\partial x^2}} \left( {{\partial u} \over {\partial x}} + {{\partial v} \over {\partial y}} \right) +
          {{\partial^2} \over {\partial y^2}} \left( {{\partial u} \over {\partial x}} + {{\partial v} \over {\partial y}} \right) \right) = LHS

Apply Continuity, so :math:`\nabla \cdot \mathbf{u} = 0`. Hence:

.. math:: -{1 \over \rho} \left( {{\partial^2 p} \over {\partial x^2}} + {{\partial^2 p} \over {\partial y^2}} \right) = LHS

The Poisson Equation in Vector Notation
---------------------------------------

Equate LHS and RHS
~~~~~~~~~~~~~~~~~~

.. math:: -{1 \over \rho} \left( {{\partial^2 p} \over {\partial x^2}} + {{\partial^2 p} \over {\partial y^2}} \right) =
           \left( {{\partial u} \over {\partial x}} \right)^2 + 
          2 {{\partial u} \over {\partial y}} {{\partial v} \over {\partial x}}+
          \left( {{\partial v} \over {\partial y}} \right)^2

In Vector Form
~~~~~~~~~~~~~~

.. math:: \nabla^2 p = -f

* A Poisson Equation for pressure, which ensures that continuity is satisfied.

* Now pressure and velocity are coupled in the continuous domain.

Poisson's Equation in Numerical Domain
======================================

We have shown that the Poisson's Equation is valid in the continuous domain, but in the numerical domain, we use discretisation

Momentum Equations in Vector Form
---------------------------------

The Momentum Equation for an **incompressible fluid**:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: {{\partial \mathbf{u}} \over {\partial t}} + \mathbf{u} \cdot \nabla \mathbf{u} =
          {-{1 \over \rho} \nabla p} + {\nu \nabla^2 \mathbf{u}} 

Discretised Momentum Equations in Vector Form
---------------------------------------------

Discretise in time: FD in time, with pressure at time :math:`n+1` (the pressure that corresponds with the velocity at :math:`n+1`)

.. math:: \mathbf{u}^{n+1} = \mathbf{u}^n + \Delta t \left(-\mathbf{u}^n \cdot \nabla \mathbf{u}^n - {1 \over \rho} \nabla p^{n+1} + \nu \nabla^2 \mathbf{u}^n \right)

The Divergence of the Momentum Equations in Vector Form
-------------------------------------------------------

Now take the divergence of the momentum equations, if :math:`\mathbf{M}` is the "vector" of momentum equations, the divergence is:

.. math:: \nabla \cdot \mathbf{u}^{n+1} = \nabla \cdot \mathbf{u}^n +
          \Delta t \left(-\nabla \cdot (\mathbf{u}^n \cdot \nabla \mathbf{u}^n) -
          {1 \over \rho} \nabla^2 p^{n+1} +
          \nu \nabla^2 (\nabla \cdot \mathbf{u}^n) \right)

* In the numerical scheme, we want a divergence free velocity at the next step, i.e. :math:`\nabla \cdot \mathbf{u}^{n+1} = 0` to satisfy continuity
* But at the current step we may have :math:`\nabla \cdot \mathbf{u}^{n} \ne 0` due to numerical errors
* So we can't cancel out the divergence of velocity in the numerical domain (although we could in the continuous domain)

This is a Fractional Step approach:

* Solve the momentum equation for velocity, numerically (the divergence of the velocity might not be zero) :math:`\mathbf{u}^{n+1/2}`
* Solve the Poisson equation for pressure, forcing the divergence to be zero
* Correct the velocity to satisfy continuity

The Poisson Equation in Vector Notation
---------------------------------------

Poisson Equation for :math:`p` at time :math:`n+1` and forcing :math:`\nabla \cdot \mathbf{u}^{n+1} = 0`

.. math:: \nabla^2 p^{n+1} = \rho {{\nabla \cdot \mathbf{u}^n} \over {\Delta t}}-
                             \rho \nabla \cdot (\mathbf{u}^n \cdot \nabla \mathbf{u}^n)+
                             \mu \nabla^2 (\nabla \cdot \mathbf{u}^n)
                             
* In the **numerical domain** the velocity field we are producing with the Navier Stokes equations is **not** completely divergence free
* Think of the velocity obtained from the Navier Stokes as being at an intermediate step :math:`\mathbf{u}^{n+1/2}`
* And :math:`\nabla \cdot \mathbf{u}^{n+1/2} \ne 0`
* We need :math:`p^{n+1}` so that **continuity is satisfied**


