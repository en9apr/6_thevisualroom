==============================
Second Order Numerical Methods
==============================

.. contents::
   :local:

Leapfrog Scheme
---------------

Equations
~~~~~~~~~

**For linear convection** the schemes we have considered so far are:

* Unconditionally unstable (1st order FD in time, **2nd order CD in space**)

.. math:: {{u_i^{n+1} - u_i^n} \over {\Delta t}} + c {{u_{i+1}^n - u_{i-1}^n} \over 2\Delta x}=0 
 
.. math:: u_i^{n+1} = u_i^{n} - {\sigma \over 2} (u_{i+1}^n - u_{i-1}^n) 
   :label: FTCS

* Introduce numerical diffusion (1st order FD in time, **1st order BD in space**)

.. math:: {{u_i^{n+1} - u_i^n} \over {\Delta t}} + c {{u_i^n - u_{i-1}^n} \over \Delta x}=0 
 
.. math:: u_i^{n+1} = u_i^{n} - \sigma (u_{i}^n - u_{i-1}^n) 
   :label: FTBS

For demonstrations of this, see numerical_scheme_1_

.. _numerical_scheme_1: http://www.thevisualroom.com/numerical_scheme_1.html

For the leapfrog scheme, both space and time are discretised by 2nd order CD formulas

.. math:: {{u_i^{n+1} - u_i^{n-1}} \over {2 \Delta t}} + c {{u_{i-1}^n - u_{i-1}^n} \over {2 \Delta x}}=0 

i.e.

.. math:: u_i^{n+1} = u_i^{n-1} - {{c \Delta t} \over \Delta x} (u_{i+1}^n - u_{i-1}^n) 

.. math:: u_i^{n+1} = u_i^{n-1} - \sigma (u_{i+1}^n - u_{i-1}^n) 

where the CFL number is :math:`\sigma =  {{c \Delta t} \over \Delta x}` 

.. figure:: ../_images/leapfrog_scheme.png
   :scale: 100%
   :align: center


* :math:`u_i^{n+1} \Rightarrow` New Solution
* :math:`u_i^{n}` Does not contribute (leapfrogging)
* 3 time levels in the discretisation
* Requires initialisation by using another method - e.g. **upwind** (a starting scheme)

von Neumann analysis
~~~~~~~~~~~~~~~~~~~~

.. math:: V^{n+1} = V^{n-1} - \sigma V^n(e^{I \phi} - e^{-I \phi})

As G is independent of n, write: 

.. math:: G = {V^{n+1} \over V^n} = {{V^{n}} \over {V^{n-1}}} 

Quadratic Equation for G:

.. math:: G - {1 \over G} = - \sigma (e^{I \phi} - e^{-I \phi}) 

Solution:

.. math:: G = I \sigma sin \phi \pm \sqrt{ 1 - \sigma^2 sin^2 \phi  }

* If :math:`\sigma > 1` the scheme is unstable, since sqrt term can be negative, thus G pure imaginary and magnitude of G > 0 (e.g. :math:`\phi = \pi / 2`)

* If :math:`\sigma < 1` them in sqrt is always real and 

.. math:: G.G^* = Re(G)^2 + Im(G)^2 = (1- \sigma^2 sin^2 \phi)+\sigma^2 sin^2 \phi = 1

Therefore the scheme is **neutrally stable** (oscillations will neither grow nor reduce) because the amplification factor is equal to 1 for the convection equation.

**Question: would this therefore be useful for analysing ocean waves that have an oscillatory input?**

Lax-Friedrichs Scheme
---------------------

Introduced in the paper:

Lax P.D. "Weak Solutions of Nonlinear Hyperbolic Equations and Their Numerical Computation", Communications on Pure and Applied Mathematics (1954)

History:

* Peter Lax laid the foundations for the modern theory of non-linear hyperbolic equations and shock wave theory. In 2005 he won the Abel Prize for mathematics. 

Equations
~~~~~~~~~

The idea of the Lax-Friedrichs scheme is to replace :math:`u_i^n` in :eq:`FTCS` by the average:

.. math:: u_i^n = {1 \over 2} (u_{i-1}^n + u_{i+1}^n) 

This will stabilize FD in t / CD in x (Forward Time, Centred Scheme - FTCS)

.. math:: u_i^{n+1} = {1 \over 2} (u_{i-1}^n + u_{i+1}^n) - {\sigma \over 2}(u_{i+1}^n - u_{i-1}^n)

Substitution introduces an error :math:`O(\Delta x) \Rightarrow` Reduces the order of the scheme to **first order** - however it is now stable (FTCS was unconditionally unstable for the convection equation)

.. figure:: ../_images/lax_friedrichs.png
   :scale: 100%
   :align: center

:math:`u_i^{n+1}` does not depend on :math:`u_i^n`

von Neumann analysis
~~~~~~~~~~~~~~~~~~~~

Insert the following into discretized equation:

.. math:: {V^n} {e^{I \phi}}

In the usual way we obtain:

.. math:: G = cos \phi - I \sigma \phi 

This results in an ellipse in the complex plane

.. figure:: ../_images/lax_friedrichs_stability.png
   :scale: 60%
   :align: center

CFL stabilty condition applies, :math:`\sigma \le 1`

**Question: is this useful for shock wave modelling, because the scheme introduces an artificial viscosity term, i.e.** :math:`1 \over 2` **?**

Lax-Wendroff Scheme
-------------------

Introduced in the paper:

Lax, P. and Wendroff, B. "System of Conservation Laws", Communications on Pure and Applied Mathematics (1960)

The Lax-Friedrichs scheme stabilized FTCS scheme, but introduced an error that was too large, i.e. **unacceptable 1st order error**. 

**The Lax-Wendroff scheme was the first scheme introduced that was 2nd order in space and time - with only TWO time levels (unlike the Leapfrog scheme which has THREE)**

**History:** This is a landmark scheme in the history of CFD and was used in aeronautical applications from the 1960s - 1980s

Equations
~~~~~~~~~

**Procedure**

* Taylor expansion in time:

.. math:: u_i^{n+1} = u_i^n + {\Delta t}(u_t)_i^n + {\Delta t^2 \over 2} (u_{tt})_i^n + O(\Delta t^3)
   :label: one

where:

.. math:: u_{t} = {{\partial u} \over {\partial t}}

* Keep the second time derivative in the discretisation
* Replace the time derivatives by equivalent space derivatives

**Application of Procedure**

* Use convection equation: :math:`u_t + cu_x = 0 \quad \Rightarrow \qquad u_t = -cu_x`

* Take the time derivative of the convection equation: :math:`\partial / {\partial t} \quad \Rightarrow \qquad u_{tt} = -c(u_x)_t = -c(u_t)_x = c^2u_{xx}`

* Replace :math:`u_t` and :math:`u_{tt}` in :eq:`one`:

.. math:: u_i^{n+1} = u_i^n - c{\Delta t}(u_x)_i^n + {{c^2 \Delta t^2} \over 2} (u_{xx})_i^n + O(\Delta t^3)

This has introduced an additional dissipative term :math:`{{c^2 \Delta t^2} \over 2} (u_{xx})_i^n`

* Using CD in x on the Taylor Expansion results in the **Lax-Wendroff Scheme**:

.. math:: u_i^{n+1} = u_i^n - {\sigma \over 2}(u_{i+1}^n - u_{i-1}^n) + {\sigma^2 \over 2} (u_{i+1}^n - 2u_i^n + u_{i-1}^n)

.. figure:: ../_images/lax_wendroff.png
   :scale: 100%
   :align: center

Notes on Lax-Wendroff Scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Looking back at the Modified Differential Equation for FTCS

* LW scheme is the discretisation of a modified convection equation obtained by adding the lowest order truncation error term:

.. math:: u_t + cu_x + {{\Delta t} \over 2} c^2 u_{xx} = 0

* LW dominating truncation error is :math:`\sim u_{xxx}` and it's modified differential equation:

.. math:: \bar{u}_t + c \bar{u}_x + {{\Delta t} \over 2} c^2 \bar{u}_{xx} = {{c \Delta x^2} \over 6} \bar{u}_{xxx} + O(\Delta t^2, \Delta x^4)


von Neumann analysis
~~~~~~~~~~~~~~~~~~~~

The result is the following amplification factor:


.. math:: G = 1 - {\sigma \over 2} (e^{I \phi} - e^{-I \phi}) + {\sigma^2 \over 2} (e^{I \phi} - 2 + e^{-I \phi})
            = 1 - I \sigma sin \phi - \sigma^2(1- cos \phi)

Real and imaginary parts:

.. math:: \xi = Re(G) = (1- \sigma^2) + \sigma^2 cos \phi
 
.. math:: \eta = Im(G) = -\sigma sin \phi

This results in:

.. figure:: ../_images/lax_wendroff_stability.png
   :scale: 50%
   :align: center

**Lax Wendroff scheme is stable if** :math:`\sigma < 1` **(CFL condition)**
