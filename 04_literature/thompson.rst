============================================================
Thompson et al. (2001) "Investigation of Rigid-Roll Coating"
============================================================

Thompson, H.M., Kapur, N., Gaskell, P.H., Summers, J.L. and Abbott, S.J. "A theoretical and experimental investigation of reservoir-fed rigid-roll coating", Chemical Engineering Science, vol 56, 2001.

.. contents::
   :local:

Abstract
========

A variant of reverse roll coating is studied via:

* Experimental
* Analytical (lubrication)
* Computational (finite element)

Experimental results:

* Flow rate and wetting line position over range of roll speed ratios and capillary numbers
* Provided the wetting line is sufficiently far from the nib, the flow rate depends linearly on reservoir level

What is unique about this study:

* Variation of dynamic contact angle with metering roll speed has been accounted for

Lubrication model boundary conditions:

* Free surface
* Surface tension
* Wetting line effects 

Physics:

* The effect of gravity is influential
* Flow in the reservoir is recirculating in nature
* The size and number of recirculations depends on reservoir geometry

Good agreement between models and experimental data

Introduction
============

Rigid-roll coating systems operating in reverse mode are used for applying thin liquid films, in the production of:

* Magnetic media
* Adhesive tapes
* Films
* Foils

Types of rigid-roll coating systems:

* Pan-fed
* Direct-fed

Physics:

* Flow in the **positive metering gap** between the co-rotating applicator and metering rolls plays a crucial role in determining the **final wet film thickness transferred to the substrate**

Discrepancy between FE model and experiment:

* Neglect of gravity
* Absence of dynamic wetting line

Discrepancy between lubrication model and experiment:

* At high speed ratio and capillary number
* Importance of the position of the dynamic wetting line relative to the nip - w.r.t onset of cascade instability

Other effects:

* Effect of feed condition on the flow structure

Direct-fed reverse roll coating is investigated. The web thickness is determined by:

* Competition between metering action of the reverse roll configuration
* Influence of gravity in the form of hydrostatic head provided by the reservoir

Aim of this paper
-----------------

Quantify the effect of the reservoir on:

* The coated film thickness
* Dynamic wetting line position


Experimental Method
===================

A rig was designed to:

* Retain the important features of the industrial model
* Be flexible in operation
* Allow optical access

Conditions:

* Large bank of liquid ensured that ambient temperature changes did not affect liquid properties
* 0.5kg weight was applied to the hosing of the tank to seat it on the rolls
* Doctor blade used to remove any residue

Fluid:

* Newtonian mineral oil
* Surface tension measured with a torsion balance
* Viscosity measured with a viscometer
* Density measured by measuring the mass of a known volume of fluid

Visualisation:

* Liquid illuminated
* Illumination of nip region achieved using a mirror
* Flow patterns recorded on a video system via a microscope and CCD camera


Flux and wetting line location measurements
-------------------------------------------

Non-contacting methods for measuring film thickness:

* Infra-red absorption
* Microwave absorption
* X-ray fluorescence
* Capacitance probes

Problem with non-contacting methods:

* Not possible to deduce an exact flux due to variations in velocity through the film

Solution:

* Web was scraped clean of liquid using a twin scraper blade
* A mass of liquid was collected over a known time interval, hence flux determined

Wetting line:

* Position of wetting line measured using microscope fitted with a cross hair graticule


Flow visualisation
------------------

Purpose:

* Highlight flow patterns near dynamic wetting line
* Show how the flow patterns depend on operating parameters

Hydrogen bubble technique impractical due to:

* buoyancy effects
* long residence times
* difficulty of constraining a bubble stream to a narrow planar section

Used dye injection method instead (a laser was not needed)

Flow in reservoir was visualised by:

* discharging dye continuously by traversing flow field at stagnation points
* dye was not confined to a single streamline

Flow in nip region was visualised by:

* releasing a pulse of tracer liquid onto metering roll just upstream of dynamic wetting line

Mathematical Models
===================

Assumptions:

* Isothermal
* 2D
* Newtonian
* Incompressible
* Governed by the Navier Stokes Equations

Two approaches:

1) Lubrication model:
  
  * A new approach - where variation of contact angle will roll speed is predicted rather than prescribed

2) Numerical model:

* Finite element method for full non-linear problem
* Algebraic mesh generation algorithm

Lubrication Analysis
--------------------

Assumptions:

* Uni-directional
* Away from free surface
* One-dimensional

Physics:

* Used non-dimensional scaling
* In the past the influence of wetting lines has been ignored in previous free surface models
* Lubrication theory is inapplicable in the vicinity of the downstream free surface and the solution domain is terminated by introducing appropriate boundary conditions

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

* At the free surface, the pressure is equal to the capillary pressure due to surface tension
* We assume that the free surface forms an arc of a circle, so it is a function of the dynamic wetting line angle
* Previous studies didn't model the idea that the dynamic wetting line angle varies with metered roll speed, however **the nature of the flow near the dynamic wetting line is still a matter of debate**
* The variation is accounted for using a hydrodynamic asymptotic model for wetting
* The asymptotic model is based on the assumption that the free surface is planar near the wetting line and requires the velocity of the fluid in the liquid-gas interface to be a function of the dynamic wetting line angle
* A calibration procedure is adopted in order to estimate the variables in the above in terms of a single adjustable parameter - the interfacial thickness. **Further experimental data is needed to verify the accuracy of the estimates for the parameters**
* The calibration results in two equations in terms of Bond Number (which measures the relative importance of gravitational to surface tension force) and Capillary Number
* A third equation is derived to relate the flux to the radius of curvature of the meniscus at the wetting line
* Newtonian iteration is used to solve the three equations

Limitation:

**Lubrication theory is inapplicable near the meniscus where the flow is 2D and often recirculating**

This requires fully non-linear 2D finite element simulation

Finite Element Solutions
------------------------

Flow in the reservoir is solved using FE method (due to it's topological flexibility)

Models:

1) Extends from top of the reservoir into the nip region
2) Extends from the nip region to the downstream free surface

Solution method:

* The governing equations are non-dimensionalised w.r.t. velocity, length and pressure
* Solved using Bubnov-Galerkin weighted residual FE formulation

Problems:

* Obtaining FE solutions in the reservoir region is straightforward since the geometry is fixed.
* Flow in the downstream region is more complex because of the existence of the meniscus and the associated wetting line

**Downstream Model**

Solution to boundary conditions for downstream model:

* Boundary conforming mesh based on spine approach 
* **Conditions at the wetting line are controversial**

Solution:

* Dynamic wetting angle is predicted and fluid is allowed to slip on the roll surface for those nodes adjacent to the wetting line

**Reservoir Model**

* Velocity profile at nip is the same as that specified in the downstream model
* Reservoir surface set to atmospheric pressure
* Parabolic velocity profile set at the surface to ensure mass conservation

**Solution**

* Frontal Method
* Newton iteration
* Leading to second order convergence

Grid refinement was also completed

Results
=======

Variables of interest:

* Flow rate
* Wetting line position is useful w.r.t cascade instability when it migrates upstream of the nip

Two experiments:

1) Peripheral speed of the applicator **fixed** and metering roll speed **varied**. Three different reservoir levels. Capillary number **fixed**, speed ratio **varied**

2) Peripheral speed of metering roll **fixed** and applicator roll speed **varied**. Capillary number **varied** and speed ratio **varied**

In both sets of experiments Bond Number is **fixed** 

1) Experiment 1

* Shows q decreasing as S increases until a critical value is reached, beyond which q increases. This is caused by the wetting line moving upstream of the nip, which increases the effective gap at the meniscus.
* Effect of reservoir level depends on position of wetting line w.r.t. nip - i.e. whether it's upstream or downstream of the nip

2) Experiment 2

* Shows q increases as Ca increases
* Wetting line is located downstream of the nip for Ca > 0.2

**Validation of model:**

* q and wetting line position are validated with good agreement  

Results also show the sensitivity of the wetting line position to dynamic wetting line angle 

**Problem** we are unable to predict wetting line positions at the higher Ca numbers - possibly due to Bretherton condition
**Suggesting** flow rate insensitive to this condition, whereas wetting line position is more sensitive.
**Reasoning** Circle assumption is not a good representation of the free surface shape for higher Ca numbers - it's probably more parabolic

Effect of head:

* linear dependence between q and head - where wetting line is safely downstream of the nip

Flow inside the reservoir:

* Recirculations can cause difficulties for the lubrication model because the rectilinear assumption is violated in such regions

Conclusions
===========

1) This is a study of the equilibrium flow in a variant of reverse roll coating, where the metering gap sits beneath a large reservoir

2) Experimental data shows:

* If wetting line is sufficiently far downstream of the nip: flow rate increases **linearly** with reservoir level
* If wetting line is close to the nip: effect of level on flow rate is **non-linear** and influences the onset of the **cascade instability**

3) Hydrodynamic model:

* Dynamic contact angle with metering roll speed is incorporated using a hydrodynamic model for wetting. 
* However, data is lacking for the hydrodynamic model, so a calibration method is proposed, in terms of the interfacial layer thickness.
* This interfacial layer thickness is calibrated by matching the prediction against theory for one data point only.

4) FE solutions can predict accurately:

* Flow rate
* Wetting line positions (even at high Capillary Numbers)

5) Lubrication model predicts the flow rate over a range of Capillary Numbers, however:

* Wetting line prediction is sensitive to the way the free surface is represented
* Free surface profile should be represented by a parabolic rather than a circular arc
* But the lubrication prediction demonstrates the benefit of incorporating the effects of free surface, surface tension and the wetting line - especially the ability to predict the flow rate minimum caused by wetting line migration through the nip

6) Experimental results vs FE solutions show good comparison

7) Flow in reservoir should not be overlooked - given the possibility of recirculations there

 


