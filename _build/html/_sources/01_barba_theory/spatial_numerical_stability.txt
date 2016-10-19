==========================================================================
Amplification of the Semi-discretised System for Space-Time Discretisation
==========================================================================

.. contents:: 
   :local:

Summary of Previous Work
========================

Recall the semi-discretised system:

.. math:: {d \mathbf{u} \over dt}  = \mathbf{S} \cdot \mathbf{u} + Q
   :label: 1

Where:

* :math:`\mathbf{S}` is the matrix representing the spatial discretisation (with BCs)
* :math:`\mathbf{u}` is the vector of nodal values
* :math:`Q` is a non-homogenous term related to the boundary conditions

It's exact solution :math:`\overline{\mathbf{u}}` has a **"Model Decomposition"**:

.. math:: \overline{\mathbf{u}}(t, \mathbf{x}) = 
          \sum_{j=1}^N \overline{\mathbf{u}_j}(t) V^j(\mathbf{x})
   :label: 2

Inserting :eq:`2` into :eq:`1`:

.. math:: {d \over dt} \left( \sum_{j=1}^N \mathbf{u}_j(t) V^j (\mathbf{x}) \right) = 
          \mathbf{S} \cdot \sum_{j=1}^N \mathbf{u}_j(t) V^j (\mathbf{x}) + Q

.. math:: \mathbf{S} \cdot V^j = \Omega_j V^j(\mathbf{x})

So the equation for the time-dependent coefficients are:

.. math:: {d \over dt} \overline{\mathbf{u}_j} = \Omega_j \overline{\mathbf{u}_j} + Q_j

The homogeneous solution was :math:`\overline{\mathbf{u}}_{jT} = c_{0k}e^{(\Omega_j t)}` (the transient)

This transient completely determines the stability of the semi-discretised system :eq:`1`

We only look at this homogeneous part - assume :math:`Q=0`

Stability Condition for Time Integration
========================================

At :math:`t=n \Delta t`, then

.. math:: \overline{\mathbf{u}_T}(n \Delta t) = \sum_{j=1}^N \overline{\mathbf{u}_{Tj}}(n \Delta t) V^j = 
                                                \sum_{j=1}^N \overline{\mathbf{u}_{j}^0} e^{\Omega_j n \Delta t} V^j


Define an amplification factor of the exact solution of :eq:`1` by:

.. math:: \overline{\mathbf{u}_T}(n \Delta t) = \overline{G}(\Omega_j) \cdot \overline{\mathbf{u}_{Tj}}((n-1) \Delta t) = 
          \overline{\mathbf{u}_{j}^0} e^{\Omega_j n \Delta t} = 
           e^{\Omega_j \Delta t} \left( \overline{\mathbf{u}_{j}^0} e^{\Omega_j (n-1) \Delta t} \right)

So the amplification factor of the exact solution is:

.. math:: \overline{G}(\Omega_j) = e^{\Omega_j \Delta t}

Stability condition:

As before: :math:`Re(\Omega_j) \le 0  \quad \forall j \quad \Leftrightarrow \quad \left| \overline{G} \right| \le 1 \quad \forall \phi \in [-\pi, \pi]`

All the properties of the time integration can be looked at separately 

We just need to look at the "Scalar Modal Equation" (dropping subscript j)

.. math:: {dw \over dt} = \Omega w

:math:`w` is an arbitrary component of the Modal Equation

Analysis Time Integration Schemes
=================================

Using  :math:`{dw \over dt} = \Omega w` is called the "Canonical Modal Equation"

Stability Regions in the Complex plane (:math:`\Omega` - plane)

Define "Time-shift operator"

.. math:: \overline{E} \Rightarrow \overline{E}w^n = w^{n+1} \Rightarrow \overline{E}^k w^n = w^{n+k} 

General Time Integration Method:

.. math:: w^{n+1} = P(\overline{E}, \Omega \Delta t) w^n

:math:`P` has the effect of being a numerical amplification factor:

.. math:: w^n = P(\overline{E}, \Omega \Delta t) w^{n-1} = \cdots = P(\overline{E}, \Omega \Delta t)^n w^{0}

where :math:`w^0` is at time=0

Stability requires that :math:`w^n` must stay bounded, i.e. :math:`P^n` uniformly bounded :math:`\forall n` and :math:`\Delta t`

In particular :math:`n \Rightarrow \infty` and :math:`\Delta t \Rightarrow 0` (with :math:`n \Delta t` fixed)

i.e. :math:`\lvert P \rvert \lt k` (independent of :math:`n` and :math:`\Delta t`) for :math:`0 \lt n \Delta t \lt T` (finite time)

:math:`z_P` are the eigenvalues of P (solutions to the characteristic polynomial), i.e. solution of

.. math:: z_P = P(z_P, \Omega \Delta t)

Neccesary condition (not always sufficient) for stability 

.. math:: \left| z_P \right| \le 1
   :label: condition

For all space discretisation that satisfy :math:`Re(\Omega_j) \le 0 \quad \forall j`

The associated numerical discretisation in time will be stable of condition :eq:`condition` is satisfied

Analysis of Space-Time Discretisation
=====================================

Compare numerical amplification factor :math:`z_P` with the exact amplification factor. We have

.. math:: w^n = [P(\overline{E}, \Omega \Delta t)]^n \cdot w^0 = z_P^n (\Omega \Delta t) \cdot w^n

:math:`z_P^n (\Omega \Delta t)` is the numerical approximation to the exact amplification factor :math:`\overline{G} = e^{\Omega \Delta t}`

Example 1 - Forward Euler
=========================

Explicit FD in time

Applied to:

.. math:: {dw \over dt} = \Omega w 

Leads to:

.. math:: w^{n+1} - w^n = \Omega \Delta t w^n

Leading to:

.. math:: z_P = P = 1+\Omega \Delta t

Therefore stable for all discretisations associated to an eigenvalue spectrum, such that:

.. math:: \left| 1 + \Omega \Delta t \right| \le 1

Or:

.. math:: [1 + Re(\Omega \Delta t)]^2 + [Im(\Omega \Delta t)]^2 \le 1

In a :math:`\Omega \Delta t` complex plane, this is a circle of unit radius centred at :math:`\Omega \Delta t = -1`:

.. figure:: ../_images/convection_spectrum_3.png
   :align: center
   :scale: 70%

Diffusion Operator
------------------

.. math:: {d u_i \over dt} = {\alpha \over \Delta x^2}(u_{i+1} - 2u_i + u_{i-1})

We previously obtained the matrix S and found the eigenvalues:

.. math:: \Omega(\phi_j) = {2 \alpha \over \Delta x^2}(cos \phi_j - 1)

i.e. the eigenvalues are real and negative

The stability condition is: 

.. math:: -2 \le -\left| Re(\Omega \Delta t) \right| \le 0

i.e.

.. math:: -2 \le -\left| {2 \alpha \over \Delta x^2}(2) \Delta t \right| \le 0

Hence: 

.. math:: 0 \le {{2 \alpha} \over {\Delta x^2}} \le {1 \over 2}

**i.e. stable**

Which is the same as we had from von Neumann (although the method is different - we have separated time and space analysis of the stability)

Convection with CD in Space
---------------------------

We previously obtained the matrix S and found the eigenvalues:

.. math:: \Omega(\phi_j) = -I {a \over \Delta x}(sin \phi_j)

So:

.. math:: \Omega \Delta t = -I {{a \Delta t} \over \Delta x}(sin \phi) = -I \sigma sin \phi

Which is purely imaginary, outside the stability circle of the Forward Euler method

**This is an unstable combination** 

Convection with Upwind
----------------------

We previously obtained the matrix S and found the eigenvalues:

.. math:: \Omega(\phi_j) = -{a \over \Delta x}(1 - cos \phi_j + I sin \phi_j)

Hence:

.. math:: \Omega \Delta t = -\sigma(1-cos \phi +I sin \phi)

In the complex plane :math:`\Omega \Delta t` is a circle centred at :math:`-\sigma` with radius :math:`\sigma`

This circle is inside the region of stability of Forward Euler, where :math:`\sigma \le 1`, i.e. we recover the CFL condition :math:`0 \le \sigma \le 1`

Example 2 - Central Time Differencing (Leapfrog Method)
=======================================================

Explicit CD in time, leads to a 3 level 2 step method

.. math:: w^{n+1} - w^n = 2 \Omega \Delta t w^n

Hence:

.. math:: P(\overline{E}, \Omega \Delta t) = \overline{E}^{-1} + 2 \Omega \Delta t

Eigenvalues :math:`z_P` from:

.. math:: z_P = {1 \over E_P} + 2 \Omega \Delta t

Implies:

.. math:: z_P^2 - 2 \Omega \Delta t z_P - 1 = 0

Two solutions:

.. math:: z_P = \Omega \Delta t \pm \sqrt{(\Omega \Delta t)^2 +1}

Recall

* Space discretisation requires eigenvalues on left hand plane for stability

* Time integration method requires :math:`\left| z_P \right| \lt 1` for all :math:`\Omega` eigenvalues for the space discretisation 

Notes: 

* Every route :math:`z_P(\Omega \Delta t)` has to remain inside unit circle
* If some roots come outside the unit circle, when :math:`\Omega \Delta t` covers its spectrum, the scheme is unstable
* A method with two or more time levels will generate more than one solution
* When this happens, the consistency of the scheme requires than one of the eigenvalues should represent an approximation to the physical behaviour - the "Principal Solution" 

How to recognise the Principal Solution: is should tend to 1 when :math:`\Omega \Delta t \Rightarrow 0`

Physical solution:

.. math:: \lim_{\Omega \Delta t \Rightarrow 0} z_{P1} (\Omega) = 1

The other solution is called the "Spurious Solution" - represents a non-physical time behaviour (introduced by the scheme)

Back to Leapfrog:

.. math:: z_P = \Omega \Delta t \pm \sqrt{(\Omega \Delta t)^2 + 1}

How do these behave as :math:`\Delta t \Rightarrow 0`:

.. math:: z_1 = \Omega \Delta t + (1 + (\Omega \Delta t)^2)^{1/2}
              = 1 + \Omega \Delta t + {(\Omega \Delta t)^2 \over 2} - {(\Omega \Delta t)^4 \over 8} + \cdots
   :label: one

.. math:: z_2 = \Omega \Delta t - (1 + (\Omega \Delta t)^2)^{1/2}
              = -1 + \Omega \Delta t - {(\Omega \Delta t)^2 \over 2} - {(\Omega \Delta t)^4 \over 8} + \cdots
   :label: two

:eq:`one` is the physical solution (because it tends to 1)

:eq:`two` is non-physical (because it tends to -1)

Recall:

.. math:: \overline{G}(\Omega) = e^{\Omega \Delta t} = 
          1 +  \Omega \Delta t - {(\Omega \Delta t)^2 \over {2!}} - {(\Omega \Delta t)^3 \over {3!}} + \cdots

Compare with :eq:`one`, the first three terms are exactly the same, so the scheme is second order in time

Characteristic Polynomial

.. math:: {z_P - {1 \over z_P}} = 2 \Omega \Delta t

With a stability limit :math:`z_P = e^{I \Theta}`

We obtain :math:`\Omega \Delta t = I sin \Theta`

.. figure:: ../_images/segment.png
   :align: center
   :scale: 70%

Conclusion

The diffusion operator and upwind convection have real negative eigenvalues

This will lead to unstable scheme when solve by Leapfrog. Leapfrog does not handle dissipative schemes.

Central differencing will be ok with Leapfrog integration

Example 3 - Backward Euler
==========================

Implicit backward difference in time:

.. math:: w^{n+1} - w^n = \Omega \Delta t w^{n+1} = \Omega \Delta t (\overline{E} w^n)

Or:

.. math:: (1+\Omega \Delta t)w^{n+1} = w^n

Hence:

.. math:: w^{n+1} = \Omega \Delta t \overline{E} w^n + w^n

Leading to:

.. math:: P(\overline{E}, \Omega \Delta t) = 1 + \Omega \Delta t \overline{E}

Eigenvalue :math:`z_P`:

.. math:: z_P = 1 + \Omega \Delta t z_P

Or:

.. math:: z_P = {1 \over {1 - \Omega \Delta t}}

Stability limit:

.. math:: z_P = e^{I \theta}

We get:

.. math:: (\Omega \Delta t) = 1 - e^{-I \theta}

Representing a circle centred on :math:`\Omega \Delta t = 1`

.. figure:: ../_images/backward.png
   :align: center
   :scale: 70%

For :math:`\left| z_P \right| \lt 1` we need :math:`\left| 1 - \Omega \Delta t \right| \gt 1`:

**ALL spatial schemes seen up to now will be stable, with implicit Euler**

Cannot look at space and time separately - can only look at space and time stability together.
