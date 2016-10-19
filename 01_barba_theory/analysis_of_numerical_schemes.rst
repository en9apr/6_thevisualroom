==================================
The Modified Differential Equation
==================================

The Truncation Error
====================

* The *truncation error* is the difference between the numerical scheme and the differential equation.

.. math::
   
   \epsilon_T =  {\Delta t \over 2} \left . {\partial^2 u \over \partial t^2} \right \vert_i^n +
          c {\Delta x^2 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i^n +
          O(\Delta t^2) + O(\Delta x^4)


* NOTE: The exact solution of the numerical scheme satisfies a MODIFIED differential equation 

Modified Differential Equation
==============================

* Consider **exact** solution of the **discretised** equation :math:`\rightarrow \bar{u}_i^n`:

.. math:: {{\bar{u}_i^{n+1} - \bar{u}_i^n} \over {\Delta t}} + {c \over {2 \Delta x}} {(\bar{u}_{i+1}^n - \bar{u}_{i-1}^n)} \equiv 0 

* From before, with :math:`u = \bar{u}`

.. math:: {{\bar{u}_i^{n+1} - \bar{u}_i^n} \over {\Delta t}} + {c \over {2 \Delta x}} {{(\bar{u}_{i+1}^n - \bar{u}_{i-1}^n)} } - 
          \left ( \left . {\partial \bar{u} \over \partial t}  \right \vert_i^n  +
          \left . {{c}} {\partial \bar{u} \over \partial x} \right \vert_i^n \right )  =
           \left .{{\Delta t} \over 2} {\partial^2 \bar{u} \over \partial t^2} \right \vert_i^n +
          O(\Delta t^2) + O(\Delta x^2)

* The above implies that for the exact solution:

.. math:: \left . {\partial \bar{u} \over \partial t}  \right \vert_i^n  +
          \left . c {\partial \bar{u} \over \partial x} \right \vert_i^n =
          -{\Delta t \over 2} \left . {\partial^2 \bar{u} \over \partial t^2} \right \vert_i^n -
          O(\Delta t^2) - O(\Delta x^2)
   :label: one

* Hence:

.. math:: \left . {\partial \bar{u} \over \partial t}  \right \vert_i^n = -
          \left . c {\partial \bar{u} \over \partial x} \right \vert_i^n - O(\Delta t) - O(\Delta x^2)

* Take :math:`\partial \over {\partial t}`:

.. math:: \left . {\partial^2 \bar{u} \over \partial t^2}  \right \vert_i^n = -
          \left . c {\partial^2 \bar{u} \over \partial x \partial t} \right \vert_i^n - O(\Delta t) - O(\Delta x^2)

* Or:

.. math:: \left . {\partial^2 \bar{u} \over \partial t^2}  \right \vert_i^n = -
          \left . c {{\partial \over \partial x}{\left (\partial \bar{u} \over \partial t \right) }} \right \vert_i^n - O(\Delta t) - O(\Delta x^2)

* i.e.

.. math:: \left . {\partial^2 \bar{u} \over \partial t^2}  \right \vert_i^n =
          \left . c^2 {{\left (\partial^2 \bar{u} \over \partial x^2 \right) }} \right \vert_i^n - O(\Delta t) - O(\Delta x^2)
   :label: two

* Substitute Equation :eq:`two` into Equation :eq:`one`. The exact solution to the numerical scheme :math:`\bar{u}` satisfies the following differential equation - called the **Modified Differential Equation**

.. math:: \left . {\partial \bar{u} \over \partial t}  \right \vert_i^n  +
          \left . c {\partial \bar{u} \over \partial x} \right \vert_i^n = -
          \left . {{c^2 \Delta t} \over 2}  {{\left (\partial^2 \bar{u} \over \partial x^2 \right) }} \right \vert_i^n -
          O(\Delta t) - O(\Delta x^2)

* Observations: The **Modified Differential Equation** is **NOT** a convection equation, it is a **convection-diffusion equation**, with a **numerical diffusion** coefficient equal to:

.. math:: {-c^2 \Delta t} \over 2

* This is negative diffusion - a process of **explosion**
* This shows why the scheme is **UNSTABLE** - it will amplify any disturbance exponentially
* The **Modified Differential Equation** and the **Truncation Error** provide essential information about the scheme

Summary of the Method for obtaining the  Modified DE
----------------------------------------------------

* Denote :math:`D(u)=0` the mathematical model we are to solve numerically (the Differential Equation)
* And :math:`N(u_i^n)=0` the numerical scheme (the Numerical Scheme)

How to obtain the Modified Differential Equation?

1) Perform consistency analysis, using Taylor series and obtain the truncation error :math:`\epsilon_T`

.. math::

   N(u_i^n)-D(u)=\epsilon_T

2) Consider the **exact** solution of the numerical scheme :math:`\bar{u}_i^n` defined by

.. math::

   N(\bar{u}_i^n)\equiv 0 

leaning to the differential equation :math:`D(\bar(u)_i^n)=-\epsilon_T`

3) Replace lowest time derivative by space derivatives in :math:`\epsilon_T` by applying Differential Equation from 2)

4) The Modified Differential Equation is defined as an equation obtained after replacement step 3), restricted to the lowest order terms (contains only space derivatives)

Example: Convection with 1st order Upwind (BD in space, FD in time)
===================================================================

* Introduce the Taylor expansions as before
* Follow steps just listed

The Modified Differential Equation looks like this:

.. math:: \left . {\partial \bar{u} \over \partial t}  \right \vert_i^n  +
          \left . c {\partial \bar{u} \over \partial x} \right \vert_i^n = 
          \left . {{c \Delta x} \over 2} \left ( 1-{{c \Delta t} \over {\Delta x}} \right ) 
          {{\partial^2 \bar{u} \over \partial x^2 }} \right \vert_i^n

The Diffusion Term is:

.. math:: {{c \Delta x} \over 2} \left ( 1-{{c \Delta t} \over {\Delta x}} \right )


CFL Condition
=============

The CFL Condition is to ensure stability of the scheme:

.. math::

  \text {For } c \gt 0

  \sigma = {{c \Delta t} \over {\Delta x}} \lt 1

:math:`\sigma` is called the CFL number, the **Courant-Friedrichs-Lewy** number
 
* CFL has a deep physical significance
* For a constant value of :math:`\sigma \lt 1` this scheme has numerical diffusion of :math:`O(\Delta x)` which is generally excessive (the scheme has poor accuracy)
