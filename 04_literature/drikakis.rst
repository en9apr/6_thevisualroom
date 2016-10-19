=========================================================
Drikakis (2003) "Advances in Turbulent Flow Computations"
=========================================================

Drikakis, D., "Advances in Turbulent Flow Computations using High Resolution Methods", Progress in Aerospace Sciences, 39 (2003) pp 405 - 424

.. contents::
   :local:

Where is this paper placed in the community?
============================================

* CFD modelling often involves attempts to reduce the effects of numerical dissipation and dispersion, but often this runs into problems with the number of grid points required and hence computing time.
* However this paper suggests that numerical dissipation and dispersion can be seen as analogous with physical dissipation and dispersion, by using the nonlinear truncation error analysis that arises from the discretisation
* This allows high resolution schemes, because instead of casting aside the numerical error, it is included in the solution.
* Although this is successful in 1D, the paper proceeds with caution about 2D or 3D solutions
* This paper seems to promote the link between physics and numerics, **because physics and numerics are intertwined** as Cavaleri et al (2007) also proposed
* It also promotes the use of numerical dissipation and dispersion depending on the conditions using limiters as Cavaleri at al (2007) suggested
* Subtle balancing of the truncation error was also something inline with Cavaleri et al (2007)

Abstract
========

* This is a review paper for high resolution methods in turbulent flow computations.
* High resolution methods successfully compute turbulent flows, without the need for an **explicit** turbulence model.
* This is a presentation of the basic properties of these methods and argues for their use as **implicit** turbulence models.
* Numerical issues that still need to be addressed: 

  - **dissipation** and **dispersion** properties, e.g. turbulence anisotropy
  - validation of under-resolved simulations of near-wall turbulent attached and separated flows

Introduction
============

What is the state of the art?
-----------------------------

* Many experimental and theoretical studies have helped to increase our physical understanding of turbulence.
* A predictive closed theory of turbulent flows has not yet been established and is unlikely to emerge in the future. 
* Hope for a **universal turbulence model** has been slowly replaced by the realization that the formulation of an adequate theory will continue to require a greatly improved understanding of the **physics of turbulent motion**.

Why is turbulence important?
----------------------------

* Lack of understanding of turbulence limits technological advancement of:

  - aircraft and car design
  - turbomachinery and combustors
 
* Also limits prediction of:

  - environmental and biological flows

What are the methods use in CFD?
--------------------------------

* Reynolds Averaged Navier Stokes Equations (RANS)
* Large Eddy Simulation (LES)
* Direct Numerical Simulation (DNS)

DNS - Direct Numerical Simulation
---------------------------------

Uses:

* Good for turbulent flow structures at **low Re** and **simple geometries**

Limitations:

* Can't resolve **unsteady flows** at **high Re**, as it's beyond foreseeable computing power
* Not sure whether DNS does provide fully resolved results - as it depends on grid convergence, which may be impractical

**RANS and LES are employed as alternatives to DNS**

RANS - Reynolds Averaged Navier Stokes Equations
------------------------------------------------

Equations are averaged over a time interval or across an ensemble of equivalent flows.

Uses:

* **Steady state** problems

Limitations:

* **Unsteadiness** introduces uncertainty.
* Requires **timescale for unsteady motion to be larger than timescale of turbulent motion** - limited to low frequency unsteady flows (which is a small number of turbulent flows).
* Closure of phase-averaged correlations is identical or very similar to conventional averaged correlations.
* Less universal than LES and DNS because of the use of experimental/DNS determination of **unknown coefficients**.
* Eddy-viscosity models can't predict **vorticity and separation**.
* However, for shock boundary layer computations, accuracy of turbulence models also depends on **discretisation of advection terms**.
 
LES - Large Eddy Simulation
---------------------------

Methodology:

* Convolve all dependent variables with a **pre-defined filter** to extract large scale components.
* All flow scales larger than filter are computed via a **modified (filtered) set of the Navier-Stokes equations (LES Equations)**.
* All flow scales smaller than filter scale (approximately the grid size) are modelled with **Sub-Grid Scale models (SGS models)**. The equations encompass unresolved correlations known as **Sub-Grid Scale stresses (SGS stresses)**.

Sources of Error in LES:

* The filter width can be constant or variable 

   - if constant additional pressure and viscous terms appear at the boundary 
   - if variable (large scale and small scale eddies) - BCs are satisfied, but LES equations now invalid - **no method exists to deal with this trade off**

* Other sources of error: 

   - numerical discretisation (dissipation and dispersion term truncation error) 
   - aliasing (in spectral schemes - arising as numerical instabilities).
   - SGS modelling (masking of SGS term by truncation error)

High Resolution Methods
-----------------------

* High resolution methods can be used to achieve the properties of Sub-Grid models, a way of **implicitly modelling turbulent flows**, known as:

   - Monotonically Integrated LES (MILES) 
   - Implicit turbulence modelling
   - Embedded turbulence modelling

* Advantages:
  
   - Characteristics implicit in these methods may mimic certain aspects of turbulence flow modelling
   - Need high resolution methods to ensure stable computation in difficult physical circumstances
 
**Physical phenomena and their numerical solution are intertwined. Question is: "how do numerical methods contribute implicitly to turbulence modelling?"** Otherwise we may double count the effect of turbulence through the explicit turbulence model as well as through the properties of the numerical method.

Investigation of embedded turbulence modelling aspects of high-resolution methods may allow understanding of CFD methods and their application in multi-dimensional problems.

Fluid Flow Equations
====================

* Mass Conservation:

.. math::  {\partial \rho \over \partial t} +
           \nabla \cdot (\rho \mathbf u) = 0
   :label: mass

* Momentum Conservation:

.. math:: {\partial (\rho \mathbf u) \over \partial t} +
          \nabla \cdot {(\rho \mathbf u \mathbf u)} = 
          -\nabla \cdot \mathbf P
   :label: momentum

* Energy Conservation:


.. math:: {\partial e \over \partial t} +
          \nabla \cdot (e \mathbf u) = 
          -\nabla \cdot (\mathbf u \cdot \mathbf P) -
          \nabla \cdot \mathbf q
   :label: energy

where:

* :math:`\mathbf u` = velocity components
* :math:`\rho` = density
* :math:`e` = total energy density per unit volume
* :math:`\mathbf q` = heat flux
* :math:`\mathbf P` = volume forces, e.g. inertial, gravitational or electromagnetic forces

The tensor :math:`\mathbf P` is defined by:

.. math::  \mathbf P = p(\rho,T) \mathbf I + {2 \over 3} \mu (\nabla \cdot \mathbf u) \mathbf I - \mu [(\nabla \mathbf u) + (\nabla \mathbf v)^T]
   :label: tensor

where:

* :math:`p(\rho, T)` = scalar pressure
* :math:`\mathbf I` = unit diagonal tensor
* :math:`T` = temperature
* :math:`\mu` = dynamic viscosity

The above system is completed by the equation of state given by:

.. math:: p = \rho R T

Irrespective of the approach employed to compute turbulent flows, the equations can be cast in the conservation laws form as follows:

.. math:: {\partial \overline{\mathbf U} \over \partial t} +
          {\partial \overline{\mathbf E} (\overline{\mathbf U}) \over \partial x} +
          {\partial \overline{\mathbf F} (\overline{\mathbf U}) \over \partial y} +
          {\partial \overline{\mathbf G} (\overline{\mathbf U}) \over \partial z} = \mathbf S(\overline{\mathbf U})
   :label: turbulence

where:

* :math:`\overline{\mathbf U}` = the array of flow variables
* :math:`\overline{\mathbf E}(\overline{\mathbf U})`,  :math:`\overline{\mathbf F}(\overline{\mathbf U})`,  :math:`\overline{\mathbf G}(\overline{\mathbf U})` = fluxes in the x, y and z directions
* The bar denotes filtering and averaging for the LES and RANS approach respectively
* DNS has no bar over the terms
* :math:`\mathbf S(\overline{\mathbf U})` = viscous terms as well as SGS and Reynolds stresses for the case of LES and RANS respectively
* For a zero RHS, we obtain a system of **hyperbolic conservation laws**
* The nonlinear terms are contained in the fluxes on the LHS

Applications of high resolution methods:

* During the last 4 decades, there have been intensive research efforts to develop accurate methods for **solving hyperbolic conservation laws**.
* The development of high-resolution methods for solving the **inviscid** equations for gas dynamics has received most of the attention.

High Resolution Methods
=======================

Properties of high-resolution methods:

* Provide at least **2nd order accuracy** in smooth areas of the flow
* Produce numerical solutions (relatively) **free from spurious oscillations**
* In the case of discontinuities, **the number of grid points is smaller** (in the transition zone containing the shock wave) in comparison with that of first-order monotone methods

Our motivation for the development of high-resolution methods emerges from our effort to circumvent Godunov's theorem, which states that:

*"There are no monotone, linear schemes for the linear advection equation of second or higher order of accuracy"*

In order words 2nd order accuracy and monotonicity are contradictory requirements.

* The key to circumvent Godunov's theorem lies on the assumption made in the theorem that the schemes are linear
* Therefore if we want to design methods which provide at least **second order accuracy** and at the same time avoid spurious oscillations in the vicinity of large gradients, then we need to develop **non-linear methods**.

* This is done in **1D** due to a lack of adequate theory for multi-dimensions.
* Even though a scheme can be designed to be 2nd order in **1D**, it might not be 2nd order in 2D or 3D.

Properties of high-resolution methods
-------------------------------------

Total Variation
~~~~~~~~~~~~~~~

Consider the one-dimensional version of Equation :eq:`turbulence` without source terms:

.. math:: {\partial \mathbf U \over \partial t} +
          {\partial {\mathbf E} ({\mathbf U}) \over \partial x} = 0
   :label: solution

Numerical approximations to weak solutions :math:`w_i` which can be obtained, for example, by :math:`(2k+1)` point explicit schemes in conservation form:

.. math:: w_i^{n+1}=w_i^n-{\Delta t \over \Delta x}(\overline{\mathbf E}_{i+1/2}^n-\overline{\mathbf E}_{i-1/2}^n)
   :label: weak_solution

where for the numerical flux :math:`\overline{\mathbf E}` yields:

.. math:: \overline{\mathbf E}_{i+1/2}^n = \overline{\mathbf E}(w_{i-k+1}^n,...,w_{i+k}^n)

and :math:`n` denotes the time level and :math:`i` is cell centered and :math:`i+1/2` is an intercell value

The numerical flux should also be consistent with the flux :math:`\mathbf E` i.e.

.. math:: \overline{\mathbf E}(\mathbf U,...,\mathbf U) = \mathbf E(\mathbf U)

Weak solution :eq:`weak_solution` converges to a solution of :eq:`solution` when the following conditions are satisfied:

1. The **total variation** of the solution w.r.t :math:`x` is uniformly bounded w.r.t :math:`t, \Delta t` and :math:`\Delta x`
2. The scheme :eq:`weak_solution` satisfies the entropy condition
3. The entropy condition implies unique solution of the initial value problem

Conditions 1. and 2. can be satisfied by the addition of **artificial viscosity** to the numerical scheme. 

* This will provide non-oscillatory solutions
* Will lose physical information (thus losing accuracy)

Total variation is defined as:

.. math:: TV(u^n)=TV(u(t))= \sum_{i=-\infty}^{+\infty} \left \vert u_{i+1}^n-u_i^n \right \vert

The function :math:`u` is assumed to be either 0 or constant as the index :math:`i` approaches infinity, in order to obtain finite total variation.

Monotonicity
~~~~~~~~~~~~

Monotonicity is defined for a scalar conservation law:

.. math:: {\partial (\rho u) \over \partial t} +
          {\partial f(u) \over \partial x} =
          {\partial u \over \partial t} +
          \alpha (u){\partial u \over \partial x} = 0
   :label: monotonicity

An important property of the weak solution of the scalar initial value problem is the *monotonicity property* according to which:

* No new local extrema in :math:`x` may be created
* The value of a local minimum increases, i.e. it is a nondecreasing function and the value of a local maximum decreases, i.e. it is a nonincreasing function.

**This the total variation** :math:`TV(u(t))` **is a decreasing function of time**

The explicit scheme of Equation :eq:`weak_solution` can be written in a shorter form as:

.. math:: w_i^{n+1}=H(w_{i-k}^n,w_{i-k+1}^n,...,w_{i+k}^n) = L \cdot w_i^n 
   :label: shorter

where:

:math:`L` is an operator. 

We say that :eq:`shorter` is *total variation nonincreasing* (TVNI) if for all :math:`w`:

.. math:: TV(L \cdot w) \le TV(w)

* The scheme in :eq:`shorter` is monotonicity preserving if the finite difference operator :math:`L` is monotonicity preserving, that is, if :math:`w` is a monotone mesh function, so is :math:`L \cdot w`.
* The scheme in :eq:`shorter` is monotone if :math:`H` is a monotone increasing function of each of it's :math:`2k+1` arguments. 

The hierarchy of these properties can be stated as follows: **the set of monotone schemes is contained in the set of TVD schemes and this is in turn contained in the set of monotonicity preserving schemes**

* Monotone schemes can be constructed as upwind or central. TVD and Essentially Non-Oscillatory (ENO) schemes can also be designed to be monotone in 1D.
* However, the set of monotone schemes is the smallest set of schemes and is a subset of TVD schemes

For a constant coefficient :math:`\alpha(u) = \alpha` we obtain the linear advection equation
Well known schemes such as the **Godunov 1st order upwind scheme** and the **Lax-Wendroff scheme** among others can cast in the general form:

.. math:: w_i^{n+1}=\sum_{l=-k_L}^{l=k_R} b_l w_{i_l}^n

where:

:math:`k_L` and :math:`k_R` are two non-negative integers and :math:`b_l` are constant coefficients. Harten has shown that the linear finite difference approximation above is monotonicity preserving if the coefficients :math:`b_l` are non-negative.

**Hence, any linear monotonicity preserving scheme is a monotone, first order accurate scheme**

Limiters
--------

What are limiters?
~~~~~~~~~~~~~~~~~~

* Limiters are the general non-linear mechanism that distinguish modern methods from classical linear schemes. 
* They are also known as **flux limiters** or **slope limiters**
* Limiters act as a nonlinear switch between more than one underlying linear method - thus **adapting the choice of numerical method based upon the behaviour of the local solution**

How to use limiters?
~~~~~~~~~~~~~~~~~~~~

* Base the analysis of the nonlinear method on the linear analysis of the available methods to be chosen by the limiter
* The limiter can also be included in the analysis providing a nonlinear truncation error analysis even for linear equations in order to achieve 2nd order accuracy.

What do limiters produce?
~~~~~~~~~~~~~~~~~~~~~~~~~

* **Limiters result in nonlinear methods even for linear equations in order to achieve 2nd order accuracy simultaneously with monotonicity**

How do limiters change the weak solution?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* 2nd order accurate 3 and 5 point nonlinear schemes of the form :eq:`shorter` can be re-written in the form:

.. math:: w_i^{n+1}=w_i^n-C_{i-1/2}\Delta w_{i-1/2}+D_{i+1/2}\Delta w_{i+1/2}
 
* *Harten's theorem* has proved that any scheme with positive values of :math:`C` and :math:`D` and with the summation of those values lying between 0 and 1 is a TVNI scheme.

What are the alternatives to Harten's Theorem?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Data compatibility condition as proposed by Roe - circumvents Godunov's theorem by constructing adaptive algorithms that would adjust themselves to the local nature of the solution
* The MUSCL approach also allows the construction of high resolution methods - nonlinear versions of these schemes that avoid spurious oscillations can be constructed by limiting the slopes in the original MUSCL scheme according to TVD constraints.

Why do we need limiters?
~~~~~~~~~~~~~~~~~~~~~~~~

* The construction of limiters does make sense not only for compressible flows that encompass discontinuities, but also for incompressible flows where large gradients occur, e.g. due to the appearance of vortices and turbulence
* Limiters play an essential role in the effective dissipation of a scheme as well as acting as a trigger for a dynamic dissipative scheme.

LES
~~~

* In dynamic LES, viscosity is adjusted locally based on whether the flow exhibits a similar structure at adjacent length scales.
* There is an implicit correspondence between some limiter forms to the dynamic SGS models in LES.
* Limiters also provide additional utility to compare local estimates of a derivative  - if these are close enough in magnitude the flow is treated as resolved - allowing the method to detect laminar flow.

Nature of algorithmic components
--------------------------------

**We must assess the similarities between high resolution methods and turbulence models for LES**

Godunov-type method
~~~~~~~~~~~~~~~~~~~

* Focus on the numerical intercell flux :math:`\tilde{\mathbf{E}}_{i+1/2}` of a Godunov-type method:

.. math:: \tilde{\mathbf{E}}_{i+1/2} = {1 \over 2}(\mathbf{E}_L + \mathbf{E}_R)-{1 \over 2} \vert \mathbf{A} \vert (\mathbf{U}_R - \mathbf{U}_L)

* This is a fairly generic and standard manner to introduce upwinding into a numerical method.
* The left and right edges can be accessed via interpolation from cell centres to edges
* This is the reconstruction step of the high resolution method

TVD method
~~~~~~~~~~

The intercell flux is given by:

.. math:: \tilde{\mathbf{E}}_{i+1/2} = {1 \over 2}(\mathbf{E}_L + \mathbf{E}_R)-{1 \over 2} \vert \mathbf{A} (1-\psi) \vert (\mathbf{U}_R - \mathbf{U}_L)

where :math:`\psi` is a limiter

An alternative formulation:

.. math:: \tilde{\mathbf{E}}_{i+1/2} = \tilde{\mathbf{E}}_{i+1/2}^{low} +
                                       \phi(\tilde{\mathbf{E}}_{i+1/2}^{high}-
                                       \tilde{\mathbf{E}}_{i+1/2}^{low})
                                       
where :math:`\phi` is a limiter

* "low" refers to first order flux (e.g. Godunov or Lax-Friedrichs)
* "high" refers to second order flux (e.g. Lax-Wendroff flux)

Flux decomposition
~~~~~~~~~~~~~~~~~~

Flux can be decomposed into terms that are:

* hyperbolic - the portion that is the sum of the local contributions (the mean flux)
* dissipative - the portion that is the difference in the variables is dissipative, with a magnitude proportional to the coefficient of numerical viscosity

TVD flux
~~~~~~~~

The TVD flux can be viewed in a similar way:

* the flux is a 2nd order centred flux, with a dissipative term than yields 1st order upwinding that is triggered by a limiter
* Hence in TVD, the action of the high resolution scheme is entirely dependent on the limiter acting on the numerical viscosity.

Non linear truncation error analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The relation of algorithmic components and flow physics can be obtained by nonlinear truncation error analysis:

.. math:: {\partial {\mathbf U} \over \partial t} +
          {\partial {\mathbf E} ({\mathbf U}) \over \partial x} = 0

.. math:: {\partial {\mathbf U} \over \partial t} +
          {\partial {\mathbf E} ({\mathbf U}) \over {\partial \mathbf U}}{{\partial \mathbf U} \over {\partial x}} = 0

Using a 1st order upwind scheme including the leading order (spatial) truncation error, we obtain the modified equation:

.. math:: {\partial \tilde{\mathbf U} \over \partial t} +
          \mathbf A (\tilde{\mathbf U}) {{\Delta \tilde{\mathbf U}} \over {\Delta x}} =
          {\Delta x \over 2} \left [
          \left \vert {\partial \hat{\mathbf E} (\hat{\mathbf U}) \over {\partial \hat{\mathbf U}}} \right \vert  {{\partial^2 \hat{\mathbf U}} \over {\partial x^2}} +
          sign \left ( {\partial \hat{\mathbf E} (\hat{\mathbf U}) \over {\partial \hat{\mathbf U}}} \right )
          {\partial^2 \hat{\mathbf E} (\hat{\mathbf U}) \over {\partial \hat{\mathbf U^2}}}
          \left ( {{\partial \hat{\mathbf U}} \over {\partial x}} \right )^2 \right ]

where the LHS of the above equation contains the terms arising from the discretisation of the differential equation

* The first term on the RHS is the 2nd order dissipation
* The second term on the RHS primarily dispersive and produces oscillations near discontinuities **like the Garden Sprinkler Effect?**

Lax-Wendroff Scheme
~~~~~~~~~~~~~~~~~~~

.. math:: {\partial \tilde{\mathbf U} \over \partial t} +
          \mathbf A (\tilde{\mathbf U}) {{\Delta \tilde{\mathbf U}} \over {\Delta x}} =
          {\Delta x^2} \left [
          - {1 \over 6}{\partial \hat{\mathbf E} (\hat{\mathbf U}) \over {\partial \hat{\mathbf U}}}{{\partial^3 \hat{\mathbf U}} \over {\partial x^3}}
          - {1 \over 2}{\partial^2 \hat{\mathbf E} (\hat{\mathbf U}) \over {\partial \hat{\mathbf U^2}}}  {{\partial \hat{\mathbf U}} \over {\partial x}}  {{\partial^2 \hat{\mathbf U}} \over {\partial x^2}}
          - {1 \over 6}{\partial^3 \hat{\mathbf E} (\hat{\mathbf U}) \over {\partial \hat{\mathbf U^3}}} \left ( {{\partial \hat{\mathbf U}} \over {\partial x}} \right )^2 \right ]

* The terms on the RHS have mixed effect although it is largely dispersive
* By using a limiter to hybridize the upwind and Lax-Wendroff method we can accentuate the dissipative effects seen in leading order spatial truncation error
* Limiters will be seen to introduce eddy viscosity
* If the limiters are not active, the dissipative term vanishes and in the case of second order methods, a fourth order dissipation results

Burgers' Equations
~~~~~~~~~~~~~~~~~~

* Margolin and Rider used truncation error analysis with numerical experiments for the Burgers' equations to show that truncation error terms for high resolution schemes have physical significance
* For non-smooth flow (e.g. in turbulence) for which the fluid velocity is not smooth over particular length and time scales, they assumed the averaged velocity is smooth at least over theses scales
* Hence the same averaged equations govern laminar and turbulent flows

Finite Difference Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Hirt made the earliest attempt to give physical meaning to numerical error - i.e. the notion that even and odd order errors are associated with diffusion and dispersion processes respectfully
* This led to the proof of entropy satisfying solutions in connection with upstream differencing.

Computational Evidence
======================

The follow examples show that high resolution methods can be used without the need to resort to SGS models.

Burgers' Turbulence
-------------------

* Burgers' equations are a 1D form of the Navier Stokes Equations
* The example used periodic boundary conditions and a random initial condition for the velocity

The comparison showed:

* Modelling the unresolved scales through a SGS model does not always improve the results
* High resolution schemes designed to satisfy the TVD condition can significantly improve the predictions without even using an SGS.
* The TVD-CB scheme without an SGS model gives the closest agreement with the DNS solution

Decaying turbulence of a homogeneous incompressible fluid
---------------------------------------------------------

* This is decaying turbulence in a triply-periodic cube
* The focus is on the nonlinearity of the convective derivatives in the momentum equation
* The simulations were carried out using the nonoscillatory forward in time (NFT) advection scheme MPDATA.

The results showed:

* Excellent agreement of the MPDATA and pseudo-spectral calculations for DNS - MPDATA is at least as accurate as the pseudo-spectral model
* Without viscous dissipation, unlimited enstrophy growth is predicted
* After the collapse of the spectral model, MPDATA is stable but inaccurate - this divergence is due to the flow topology condition being enforced by the MPDATA scheme (e.g. an embedded SGS model develops) - enstrophy has increased to the point where velocity gradients are so large that the local derivatives should be limited to obtain stability of the computations - **MPDATA results in an effective viscosity**
* High resolution schemes may be considered as producing the most stable LES result for an inviscid flow at a given computational resolution

Convective planetary boundary layer
-----------------------------------

* Previous studies used an atmospheric code based on MPDATA to accurately produce (in agreement with field/lab data and benchmark computations) the structure of the **convective planetary boundary layer**
* When an explicit turbulence model was implemented, MPDATA did not add unnecessary diffusion.
* When no explicit turbulence model was employed the high resolution model appeared to contain an effective SGS model.
* With an explicit turbulence model and when the eddy viscosity is reduced by some factor, MPDATA adds just enough dissipation.
* **These numerical experiments show the self-adaptive nature of the high-resolution methods and suggest the physically realistic character of it's truncation error**

Validation with LES with and without SGS

* The accuracy of the high resolution method in LES means a SGS model isn't needed.
* In contrast to linear methods, the success of high resolution schemes is due to the **self-adaptiveness** of the scheme during the simulation:

  - with SGS included, the resolved flow is smooth and high resolution is off
  - with no SGS, high-resolution scheme adapts the numerics assuring solutions that are as smooth as SGS

**The dissipation of high-resolution methods cannot be universally quantified since the advection scheme can be non-dissipative or dissipative, depending on the presence or absence of an explicit SGS model**

* **Open Question:** "The limits of the embedded turbulence modelling approach in the context of under-resolved simulations for wall-bounded flows need to be further investigated"

* **Open Question:** "In coursely resolved simulations, the presence of unresolved boundary layers in the flow may require explicit models to account for wall forcing"

Compressible Open Cavity Flow
-----------------------------

* The results in the preceding sections suggest that nothing really prevents implementation of high resolution methods for flows with wall boundaries (**perhaps apart from the need for finer grid resolution in the vicinity of those boundaries**).
* Use of **high resolution methods** can offer better stability and accuracy in coarsely resolved flows even with no turbulence model.
* Hybrid schemes have been used for open cavity flows - combining a Riemann solver and a flux splitting scheme

**the most dissipative scheme carries more weight in the vicinity of discontinuities where more dissipation is required**

At present in industry, under-resolved flows are used because of short turnaround times. In order to compute the compressible turbulent flow, we need to capture:

* large and small vortical structures
* free shear layers
* transitional flow
* flow separation
* flow re-laminarisation
* shock and rarefaction waves

The induced oscillatory pressure field from the cavity depends on:

* geometry
* inflow boundary layer
* freestream velocity

Previous work has used RANS (:math:`k-\epsilon` models) and DNS, but the averaging in RANS may cause concern because **turbulent timescales are much smaller than convective timescales**.

Cavity flows are characterised by:

* unsteady boundary layer separation
* instabilities within the shear layer
* dominant acoustic flow oscillations

An increase in the freestream Mach number also induces "switching" of modes, e.g. shear layer mode and wake mode (and is a function of Mach number)

Shear layer mode occurs at the downstream cavity edge and is subjugated by a single frequency.

When no explicit SGS is added, a high resolution method will not add unnecessary diffusion.

Shock-bubble Interaction
------------------------

* In the case of the inviscid solution the shrinking/disappearance of the bubble material is due entirely to numerical diffusion - means it's difficult to compare experiment and simulation
* The details of the simulated flow depending on the numerical scheme employed
* There are strong similarities between the solutions of certain methods
* Some schemes are less diffuse than others

Relation of numerics with the physics of turbulent flow and it's models
=======================================================================

Kolmogorov spectrum
-------------------

* The **Kolmogorov spectrum** describes how the **energy density** of turbulent structures decreases rapidly with increasing wave number.
* The **Kolmogorov scale** is the scale at which the viscous dissipation dominates the inertial flow
* The **turbulent cascade process** is the downward transfer of energy from large to small scales. This stops at the Kolmogorov scale, where an eddy is so small it diffuses rapidly.

* Kolmogorov defined a dissipation of kinetic energy that was independent of the coefficient of viscosity in the limit of infinite Re. 
* The **average time rate of change of dissipation of kinetic energy**, :math:`K` is given as:

.. math:: \langle K_t \rangle l = {5 \over 4} \langle (\Delta u)^3 \rangle
   :label: average


* In homogeneous, isotropic turbulence, this term is proportional to the average velocity difference at a length scale, :math:`l` cubed. 
* Note that this theory is analytic and independent of viscosity 
* This theory provides a basis for the functional form of nonlinear eddy viscosity.

Energy transfer
---------------

* The physics of the **turbulent cascade** is controlled by **the process of dissipation of this energy due to molecular viscosity at macroscopic scales** (larger than the Kolmogorov scales)
* Energy transfer is by **local eddy interactions** - energy extraction happens between energy levels of an order of magnitude difference

Shocks
------

Bethe derived the dissipation rate due to the passage of a shock:

.. math:: T \Delta S = {{G \rho^3 c^2} \over 6} (\Delta V)^3
   :label: dissipation

For Burgers' equation, a similar result may be obtained:

.. math:: \langle K_t \rangle \mathcal{L} = {1 \over 12} \langle (\Delta u)^3 \rangle 

* This is an entropy condition for Burgers turbulence, describing the minimum integral amount of inviscid dissipation for a physically meaningful solution.
* This dissipation is produced at the shocks and is a consequence of and proportional to the jump in dependent variables.

Eyink says that the dissipation of kinetic energy as defined by the Kolmogorov similarity is both **local and integral** - numerical methods can produce:

* a finite rate of dissipation independent of viscosity
* a local rate of dissipation

Use of high resolution methods to reproduce energy dissipation
--------------------------------------------------------------

* High resolution methods for hyperbolic conservation laws have been developed to capture with high accuracy the variation of flow and thermodynamic variables occurring across:

  - shock waves
  - other discontinuities

* These methods aim to correctly reproduce **shock dissipation** :eq:`dissipation`
* On the basis of the similarity of Equation :eq:`dissipation` and :eq:`average` it could be argued that high resolution methods also capture **Kolmogorov dissipation** :eq:`average`

Use of the entropy condition to decide physical solutions
---------------------------------------------------------

* **Lax-Wendroff scheme** ensures exact local conservation with proper entropy production
* This states that if there is a method in discrete conservation form that converges, it converges to a **weak solution**
* The **entropy condition** is necessary to choose the physically meaningful weak solution among the infinite set of possible weak solutions

Similarity between artificial viscosity and eddy viscosity
----------------------------------------------------------

* Eddy viscosity from Smagorinsky simplifies to von Neumann-Richtmyer artificial viscosity in one dimension
* The gradient of the von Neumann-Richtmyer artificial viscosity is:

.. math:: {{\partial \sigma} \over {\partial x}} =
          {\partial \over {\partial x}} 
          \left [
          -c(\Delta x)^2 {{\partial u} \over {\partial x}} 
          \left \vert {{\partial u} \over {\partial x}}  \right \vert      
          \right ]

The Smagorinsky model for the SGS "stress" in one dimension leads to:


.. math:: {{\partial \sigma} \over {\partial x}} =
          {\partial \over {\partial x}} 
          \left [
          -c_s {{\partial u} \over {\partial x}} 
          \left \vert {{\partial u} \over {\partial x}}  \right \vert      
          \right ]

**However, these terms differ in multiple dimensions**

* von Neumann viscosity can be extended to multiple dimensions, but the resulting non-linear viscosity is generally **anisotropic**, but isotropic models for artificial viscosity do exist.
* Other developments turn viscosity off should the flow field be considered to be resolved by numerical derivatives.
* There are also similarities between the Camassa-Holm equations (for dissipation) and entropy production

Similarity between LES subgrid models and high resolution models
----------------------------------------------------------------

* High resolution models have an effective subgrid model that is inherently local
* Algebraic form of high resolution models is similar to LES subgrid models with nonlinear eddy viscosity

Open Questions
--------------

* The theory for hyperbolic conservation laws in 1D is quite well developed
* Open questions exist in 2D and 3D - e.g. are some well posed, i.e. stable?
* 2D solutions involving a vortex sheet show a progression of greater and greater complexity as the mesh is refined.

Conclusions
===========

* High resolution methods are **implicit** and are motivated by the fact that almost all practical computations in engineering are **under-resolved**. 
* Numerics also play a role in terms of accuracy and efficiency - e.g. **numerical dissipation** regularises the flow, allowing shock propagation to be captured realistically, even if not fully.

Two competing criteria
----------------------

* High accuracy
* Protection against catastrophic failure due to nonlinear wave steepening or unresolved features

Methods:

* Nonlinear mechanisms (limiters) in high resolution methods guard from such catastrophic failures by triggering entropy production.

Two key questions
-----------------

* What criteria should be used to design the nonlinear mechanism that triggers entropy production?
* To what extent numerical dissipation accounts for turbulent flow effects? **Ideally we would like to quantify the numerical dissipation that is added to the computations**

From 1D to multi-dimensional
----------------------------

* The theory of numerical methods for hyperbolic conservation laws has made significant progress in 1D.
* However, we need to further understand the nonlinear behaviour of numerical methods in multi-dimensional problems
* This non-linear behaviour is closely related to the numerical mechanisms underlying the formation of **spurious solutions** in under-resolved flows.

Why do we get spurious solutions?
---------------------------------

* Generation of spurious vorticities may depend solely on the advection scheme - how numerical dissipation is partitioned between terms in the advection scheme. **However, the exact numerical mechanism is not understood. Why certain schemes show spurious solutions and others do not is an open question**
* The definition of the advective velocities in the equations can introduce a **truncation error vorticity source**. However, vorticity analysis of nonlinear approximations is difficult.
* Spurious solutions in under-resolved flows also depend on a **delicate balance of truncation errors** due to wave-speed dependent terms. This balance needs to be understood.

High Reynolds Number Near Wall Flows
------------------------------------

* Further validation for near wall flows is required for high resolution methods
* LES based on an explicit turbulence model poses challenges in high Re near wall flows, especially in the case of separation from a gently curved surface. 

