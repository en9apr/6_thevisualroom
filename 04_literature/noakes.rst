======================================================
Noakes et al. (2002) "Streak-line Defect Minimization"
======================================================

Noakes, C.J., Gaskell P.H., Thompson, H.M. and Ikin J.B. "Streak-line defect minimization in multi-layer slide coating systems", Trans IChemE Part A, vol 80, 2002.

.. contents::
   :local:

Abstract
========

This paper looks at the computational and experimental minimisation of **streak-lines** in **slide coating** systems.

The Problem: Eddies occur in critical regions which cause streaklines

The Solution: Minor geometrical adjustments and/or changes to wetting characteristics - especially in relation to static contact angles can:

* eliminate the occurrence of eddies
* reduce significantly associated surface deposit growth
* avoid formation of streak lines

Introduction
============

Problem for coating - avoid defects in:

* Liquid preparation
* Deposition
* Drying

The source of the defect:

* The bulk flow
* Instability - forced (uneven pumping) or self-excited (ribbing)

Stable operating conditions :math:`\Rightarrow` deposit-generated defects

Liquid layers are formed by:

1) Pumping liquid out of two feed slots down an inclined surface
2) Liquid being transferred from the slide lip to a moving web/substrate, stabilized by a small suction pressure 

Avoid:

* Streak line defects
* Band defects

Streak-lines caused by:

* Nicks/non-uniformities in the coating head
* Particles trapped in feed slots
* Build up of deposits 

Banding caused by:

* Uneven flow
* Poor web handling
* Presence of a wavy contact line

Fluid dynamics reasons for streak-lines:

* Existence of eddies 

Eddies cause:

* Trapped foreign particles
* Bubbles
* Solid deposits

Eddies occur in:

* Bead forming region
* In static contact line

Eddies affected by:

* Slide inclination
* Magnitude of static contact angle
* Cleaning of the slide lip
* Lip's radius of curvature

Particle trapping caused by:

* Lack of interfacial tension
* Break off of solid deposits at upper exit slot

Correlation exists between:

* Streak-line defects
* Surface deposits
* Flow structure in vicinity of exit slots

Aims of current study
~~~~~~~~~~~~~~~~~~~~~

1) Link between deposit-generated streak-line defects and nature of associated underlying flow structure
2) Provide data to assist in predicting whether defects will occur

Physics Investigated:

* Single and multi-layer flows
* Effect of slot exit
* Slide lip geometry
* Phenomena of back-wetting caused by static contact line migration past the upper slot corner

Experimental Methods
====================

* Single layer slot flow :math:`\Rightarrow` upper slot exit

Experimental results:

1) Free surface profiles
2) Visualization of streamlines
3) Velocity distributions

Two dimensional flow apparatus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Geometry variation:

* Square edge
* Chamfered edge

Method for visualization:

* **Hydrogen bubble technique**
* Local velocities measured by including a chopper in the illumination
* Deposit growth mimics by injecting a hardening agent into the flow

Cylindrical flow apparatus
~~~~~~~~~~~~~~~~~~~~~~~~~~

Variables:

* Surface finish
* Top corner profile
* Step height

Method for visualising specific regions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Bouncing probe method** used to calculate surface profile at slot exit
* Upper meniscus: Same method used in bead forming region between lip and moving web
* Lower meniscus: Illumination with a sheet of light from a laser diode source
* Profile of interface was determined by casting a shadow of a surface bubble track (the lower layer scattered light, which helped)
* Streakline formation was assessed by: coating many rolls under known conditions and recording the average number of sharp streak-lines per unit width.

Mathematical Modelling
======================

Overall equations and assumptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assumptions:

* Steady
* Isothermal
* 2D
* Incompressible
* Inelastic
* Coating liquids are often **shear thinning** however these effects are small - so assuming a Newtonian fluid is valid

Equations:

* 2D Navier-Stokes
* Continuity Equation

Dimensionless Numbers:

* Reynolds Number
* Stokes Number

Complications:

* Presence of free surfaces
* Internal interfaces

**For multi-layer systems, only a double layer system is considered as it shows the key features, namely:**

* Upper and lower free surfaces
* One internal interface
* Presence of both static and dynamic contact lines

**Approach is a Galerkin weighted residual formulation** 


Single-layer slot exit flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Geometry:

* Chamfered downstream corner
* Raised back slot
* Modified slot exit channel

Non-dimensional version of Navier-Stokes, length scale and velocity scale are:

* Slot width
* Flow rate per unit width / Slot width

The velocity profile at the inlet is given by:

* **1D form of Navier-Stokes Equations**

Outlet boundary obtained using:

* Stress evaluation method

Two cases:

1) Static contact **line** is pinned at the corner, while static contact **angle** at slot exit is predicted
2) Static contact **angle** is specified, while static contact **line** is predicted

Double-layer slot exit flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~

More complex than single-layer slot exit flow, it has:

* Internal interface - the nature of the boundary conditions which apply at an internal static contact line still remains a matter of debate
* Inter-layer diffusion occurs

Assume:

* Diffusion process is negligible
* No interfacial tension - velocity is continuous, but pressure discontinuous

Difficulty:

* Lack of interfacial tension results in an **over-specified problem** if the contact angle is imposed

Solution 1:

* If interfacial contact line is known to locate **close to the corner**, the preferred course of action is to pin it there
* If interfacial contact line locates away from the corner, solutions include:

   - geometric extrapolation
   - imposition of pressure continuity
   - imposition of **small interfacial tension** which decays to zero within a short distance downstream, this enabling the specification of a **contact angle**

Macroscopic flow field is independent of the choice of method, however **adding small interfacial tension was chosen**

Boundary conditions:

* Velocity profile at slot inlet is obtained as for single layer case
* Outflow determined by using stress evaluation method for 2 layers
* Inflow at upper slot determined from solution of 1D Navier Stokes solution

Bead forming flow
~~~~~~~~~~~~~~~~~

Boundary conditions:

* Inlet velocity profile obtained from 1D solution to Navier-Stokes equations
* Zero traction condition applied at outlet

Physics:

* Simplification is possible via lack of internal contact line
* Lower free surface meets slide feed device at static contact line
* Lower free surface meets moving web at dynamic contact line
* Static contact line is rarely pinned to the tip of the slide feed device, so is computed
* Static contact angle is prescribed

Difficulty:

* No slip hypothesis breaks down near dynamic contact angle

Solution:

* Dynamic contact angle is prescribed
* This allows slip between the web and the liquid close to the dynamic contact line

Results
=======

Single layer slot exit flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eddy physics:

* The top, layer forming slot can be responsible for disturbances
* Flow near the static contact line would be recirculating in nature
* Agreement between theory and experiment is good
* A slow recirculation near the contact line is predicted when a high backstep is used

Experimental proof of eddy presence causing streak lines:

* Gelatin containing hardening agent introduced into flow
* Hardened deposit found in eddy generating region
* Mechanism for generating streaklines was the gradual formation of brittle deposit due to long residence time associated with eddy presence.
* Deposits were dislodged by periodic cleaning to form minute fragments, which then produced streaklines

Back wetting: **migration of contact line upstream**:

* Even a small increase in flowrate causes back wetting

Implications of back wetting:

* Can cause the contact line to be other than straight
* Can cause recirculation near the static contact line
* Both of the above are potential for defect formation

Axi-symmetric rig:

* **Question** What are the maximum sustainable flow rates before back wetting occurs for three different slot exit designs?

**Method 1** in axi-symmetric rig:

* Cut back the top face to form a sharper edge
* Use a hydrophobic coating on the top face

Advantages of this:

* Approximately double the maximum sustainable flow rate before back wetting occurs

Disadvantages:

* Sharp top edge is impractical, due to hazard and difficult machining

**Method 2** in axi-symmetric rig:

* Height of slot is reduced
* Upstream of slide coated to prevent back wetting

Advantages of this:

* Causes static contact angle to increase to about 90 degrees
* Substantial reduction in deposit growth

* **Question** How is the flow structure affected by flow properties?

* Reducing viscosity increases eddy size near downstream slot exit
* Increasing flux also increases eddy size

Operability diagram:

* Parameters: Reynolds Number and Capillary Number
* Indicates operability window for a particular slot geometry
* Desirable range: Low Re and Medium to High Ca

Operability diagram :math:`\Rightarrow` to coat a low viscosity liquid, the slot geometry must be changed to enable eddy free regime

Geometric changes:

1) Chamfered slot exit :math:`\Rightarrow` recirculations are going to be small
2) Stepped slot:
   
   * Decrease in liquid velocity and its associated momentum due to widening the slot in turn leads to a reduction in size of eddies

   * Widening slot on upstream side is better than widening slot on downstream side

   * Sharp step may generate eddies :math:`Rightarrow` widening slot gradually or using chamfer is better

Slot exit flow for two merging layers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Comparison with numerical simulations and:

1) Two layer flow experiments for the case of a chamfered slot with identical layers of aqueous glycerine

2) Free surface and internal interface profiles for flow out of a chamfered slot

Good comparison.

A parametric study showed that fluxes and viscosities have greatest influence on flow structure:

* When upper flux :math:`\lt` lower flux, then increases in upper flux cause minimal change
* When upper flux :math:`\lt` lower flux, then eddies occur upstream

Geometric study:

* Small chamfer reduces size of downstream eddy
* Larger chamfers reduce size of upstream and downstream eddies to undetectable sizes

Curved diffuser:

* Flow is free from eddies
* Layer thicknesses change smoothly without humps observed previously

**Why don't they use curved diffusers? Too expensive? Hard to clean/manufacture?**

**Maybe it's the sharp edge - could be a hazard?**

Bead forming flow
~~~~~~~~~~~~~~~~~

Validation of numerical model against:

* Experimental upper and lower free surface profiles

Eddies can occur in bead forming region and are increased by:

* Coating gap increase
* Suction applied to lower meniscus is increased

Eddies cause:

* Prolonged residence time of coating liquid within the recirculation region, resulting in deposits forming at the lip face
* The deposits become grazed and break off by the passage of the web splices
* This leads to sharp streak lines

The model requires:

* Sufficiently refined computational mesh in the lower free surface region

Physics:

* Correlation between primary eddy strength and number of sharp streaklines per metre width

Affected by:

* Static contact angle
* Lip radius

Practical advice:

* Maintain a lip free of mechanical imperfections
* Ensure the surface properties of the lip are such that a high static contact angle is possible
* Lip radius is also an important factor

Conclusions
===========

1) Link between:

* surface deposit growth at both slot exits and the bead forming region 

and..

* streak line defects

2) Eddies close to the static contact lines are the cause of the deposit growth

3) Two key parameters that determine flow structure:

* Liquid viscosity
* Flux

In slot exit region
~~~~~~~~~~~~~~~~~~~

* To minimise/prevent eddies modify slot exit geometry using curvature, chamfers and/or a diffuser

In bead forming region
~~~~~~~~~~~~~~~~~~~~~~

* An area close to static contact line at the junction of the lower free surface and underside of lip surface is where eddies form
* Eddies will exist near the static contact line if the static contact angle is less than 40 degrees :math:`\Rightarrow` regular cleaning of the surface of the underside lip region is needed to avoid solid deposition
* Or increase static contact angle with a pre-surface treatment, thus weakening eddy structure.


