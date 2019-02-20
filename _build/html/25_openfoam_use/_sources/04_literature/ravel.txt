=====================================================================
Ravel et al. (2009) "Numerical Simulation of Viscous Nonlinear Waves"
=====================================================================

Ravel, A., Wen, X., Smith, M.H., "Numerical Simulation of Viscous, Nonlinear and Progressive Water Waves", Jounral of Fluid Mechanics, 637 (2009) pp 443 - 473

.. contents::
   :local:

Where is this study placed within the wave modelling community?
===============================================================

* This study is in the category of models attempting to solve the problem of non-linear wave effects in shallow water. 
* It takes the approach of deterministic (phase resolving) wave modelling in the time domain.
* It uses one phase per realisation, and seven realisations are performed (which may have been limited by the large computing time needed) but at different depths and frequencies.
* Unlike other models in this category, it also includes viscous effects, which allows energy decay rate to be observed, which is one of the open questions in this field.
* It also postulates various mechanisms for this decay rate, to try to explain the physics behind the energy decay rate difference between intermediate and deep water cases.

Summary of Findings
===================

Near the water surface:

1) Clockwise and anticlockwise rotation of the fluid at the trough and crest
2) Thicker vorticity layers and larger magnitude of vorticity compared with:

 * Inviscid, rotational flow (the Gerstner wave)
 * Low Re viscous flow (from Kinsman (1965) and Lamb (1932))

3) Thick shear layer with larger shear stress (-ve and +ve shear stress below crest and troughs of the wave)
4) Good comparison between wave energy decay rates (numerical, experimental and theoretical)

Literature Review
=================

3 categories of techniques:

1) Potential flow - inviscid, irrotational, no turbulence (analytical solution via Laplace equation for velocity potential plus boundary conditions):

   * Inaccurate for flows with a boundary layer (e.g. atmospheric boundary layer, bottom boundary layer)
   * Vorticity is zero, so fluid loses it's infinite degrees of freedom (so the fluid ceases to be a fluid in reality)

2) Inviscid rotational flow:

   * Gersther's trochodial wave theory
   * Boundary integral method - periodic and solitary waves with constant vorticity - limited due to special application of vorticity distribution
   * Boundary integral method - interaction of small amplitude waves with bottom ripples

3) Viscous flow:

   * Kinsman and Lamb (low Re): No convection terms, which linearises the Navier Stokes Equations, such that an analytical solution can be found
   * Behroozi: conservation of energy for fluid viscosity and decay relation
   * Wang & Joseph: viscous potential flow and viscous correction potential flow for decay of free gravity wave due to viscosity
   * Dutykh & Dias: Visco-potential free surface flow
 
What kind of numerical methods are there?
-----------------------------------------

1) Front tracking method
2) Boundary integral method
3) Phase field method
4) Second gradient method
5) Level set method
6) Marker and cell method
7) Smooth particle hydrodynamics
8) Finite analytic method and modified marker and cell method, using zero shear stress at the surface for velocity and energy decay rate
9) Kinematic boundary condition at the surface - zero shear stress, zero pressure or total zero shear stress on water surface - no air considered only water
10) Volume of fluid method - Hirt and Nichols (1981) for shape of free surface. Zwart, Burns and Galpin (2007) - VOF - define volume fraction of air and water separately with high resolution scheme for mass conservation in air and water, v, p and aw, af solved simultaneously in coupled system. Advantage of VOF = includes air and water to retain mass continuity and eliminates need to reconstruct the shape of the water surface during the calculation, like MAC, SPH and height function methods

What phenomena are we testing for?
----------------------------------

* Waves are transient, non-linear, rotational and viscous and often simplified to linear, irrotational and inviscid - this removes effect of vorticity and shear stress, which are strong indicators of rotational and viscous behaviour.
* Past vorticity work - velocity as a function of vorticity via Biot-Savart integral, eliminating pressure from the formulation. Zero shear stress - no wind at surface.

Do we need vorticity?
---------------------

* Yes - parasitic capillaries cause it - orbital vorticity below waves with wind - causes Langmuir circulations or turbulence balance - strong consequences. Vorticity is proportional to gradient of velocity - applicable for high Re flows. Action of viscous force shown in high vorticity and shear stress.

Questions arising from literature review
----------------------------------------

1) Is there vorticity and shear stress in a progressive wave without wind? (air following water)
2) If so, what is the maximum value of vorticity and shear stress in the water? Are they symmetric under crest and trough?
3) Is there an effect of water depth on vorticity and shear stress?

Method
======

Assumptions
-----------

1) 2D
2) Viscous fluid
3) Non-linear
4) Non-breaking
5) Intermediate and deep water - with a slope
6) Zero wind velocity - air follows wind
7) **Laminar flow**

Inputs
------

* :math:`L` = wavelength
* :math:`T` = periodic time
* :math:`a` = wave amplitude
* :math:`k = 2 \pi / L` = wavenumber
* :math:`\sigma = 2 \pi / T` = angular frequency
* :math:`c = L / T = \sigma  / k` = angular frequency
* :math:`2a / L` = wave steepness
* :math:`h` = depth in water
* :math:`h'` = depth in air

Outputs
-------

1) Velocity field
2) Streamlines
3) Vorticity
4) Shear stress

Physics
-------

**Boundary Conditions are:**

* Walls top and bottom
* Opening at inlet and outlet - non-linear, inviscid flow solution at inlet
* A 1:15 slope to dissipate the wave energy on the beech

**2D Equations:**

* Continuity equation for water
* Continuity equation for air
* Momentum equation for mixture in x
* Momentum equation for mixture in y
* Pressure-velocity coupling in x
* Pressure-velocity coupling in y

**Unknowns:**

* Volume fraction of water (scalar) (x,y) - by volume fraction constraint, we get the volume fraction of air (x,y)
* Velocity of water (x)
* Velocity of water (y)
* Velocity of air (x)
* Velocity of air (y)
* Pressure of mixture (scalar) (x,y)

**Initial Conditions:**

* Non-linear potential flow solution

**Correction for outward flow:**

* Asymmetry of flow means more flow moves under wave crest than out under trough
* Hence, slope is moved away at a velocity found by integrating the potential flow solution

Numerical Discretisation
------------------------

Continuity Equation:

* Second order backward Euler for transient term
* First/second order blended scheme for advection term

Momentum Equation:

* Second order backward Euler for transient term
* Second order scheme for convection term, diffusion term and pressure gradient term

Pressure-velocity coupling:

* Rhie-Chow interpolation

Grid Discretisation
-------------------

* Structured mesh
* 16 grid points over wave height
* 100 grid points over wavelength

Solver
------
 
* Algebraic multi-grid method

Test Cases
----------

All these waves are:

* Non-linear, i.e. :math:`ka > 0.1`
* Non-breaking, i.e. :math:`2a/L \le 0.08`

======== ================= ================== ============ ===================================================
Case     2a/L (steepness)  h/L (depth ratio)  T (period)   Description
======== ================= ================== ============ ===================================================
1 (C1)   0.04              0.2                0.7592       not steep, intermediate depth 
2 (C2)   0.04              0.6                0.6          not steep, deep
3 (C3)   0.06              0.2                0.7592       steep, intermediate depth
4 (C4)   0.06              0.6                0.6          steep, deep
5 (C5)   0.08              0.2                0.7592       very steep, intermediate depth
6 (C6)   0.08              0.6                0.6          very steep, deep
7 (CEX)  0.06              0.44               0.7          for comparison with experiment and energy density
======== ================= ================== ============ ===================================================

Results
=======

Energy Decay Rate
-----------------

* Determine Kinetic, Potential and Total Energy Decay Rate
* Compute decay rate :math:`\beta` by fitting this equation:

.. math::

   {\overline{E} \over \overline{E}_0} = \overline{\alpha}_T = \exp(-\beta x)

* Decay rate is higher for deep water case than for intermediate water case

**The rest of the paper discusses why by deep water has a higher energy decay rate by attempting to describe the physics** 

Free Surface Profiles
---------------------

Comparing the free surface profiles of:

* Deep water case
* Intermediate depth water case
* Nonlinear Stokes wave

Shows:

* Deep water case is **Symmetric**
* Intermediate depth water case is **Asymmetric**

Hence there is an effect of the bottom surface, forcing water to move towards the crest in intermediate depths

**Why does the fluid move upwards towards the crest in the case of intermediate depth?**

Velocity field
--------------

Intermediate water case shows higher peak velocity in the crest than the deep water case.

Streamlines
-----------

Fluid moves to where it encounters least resistance, this point is higher for intermediate water than deep water.

Vorticity
---------

* Vorticity distribution is oscillatory with depth
* The oscillations are caused by viscosity
* The maximum vorticity is larger and has a thicker layer for intermediate compared with deep cases
* However, the deep depth has more oscillations in the vorticity field, because it is deeper

Shear Stress
------------

* Intermediate depth case also has larger shear stress than deep case, because the bottom boundary restricts motion of the water, providing larger velocity gradient, hence larger shear stress.
* The deep water case has a bottom boundary located deeper, so the fluid at the surface is less influenced by viscous effects

Effect of Wave Steepness
------------------------

* Wave steepness also increases vorticity and shear stress because the velocity gradients are higher.

**The effect of the oscillatory vorticity field with depth seems to be a fundamental phenomenon and affects the shear stress and velocity field**

**The interface between the air and water is the most critical region in the flow, due to the high viscous shear stress created there**

