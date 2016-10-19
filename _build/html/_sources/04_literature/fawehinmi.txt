=========================================================================
Fawehinmi et al. (2005) "Experimental and CFD Analysis of Drop Formation"
=========================================================================

Fawehinmi, O.B., Gaskell, P.H., Jimack, P.K., Kapur, N., Thompson, H.M. "A combined experimental and computational fluid dynamics analysis of the dynamics of drop formation", Proc. IMechE Part C; J. Mechanical Engineering Science, vol 219, 2005.

.. contents::
   :local:

Abstract
========

This paper is a comparison between CFX and FLOW3D for the modelling of the dynamics of **drop formation** for:

* **Primary drop volume as a function of flowrate** (which is of significant engineering interest)
* **Interfacial features** (which seems less important in engineering, but may have an influence in extreme cases)

Primary drop volume
-------------------

FLOW3D is more efficient than CFX because FLOW3D uses a single phase approach, whereas CFX uses two phases. Both are able to predict the primary drop volume.

Fine interfacial features
-------------------------

For low flowrates and viscosity there are two phenomena:

* Fine interfacial threads
* Satellite drops

Both CFX and FLOW3D are poorer at predicting these, because of:

* **Free surface smearing** introduced by the VOF method.
* **Surface tension** becoming dominant

Literature Review
=================

Applications for the dynamics of drop formation
-----------------------------------------------

* Separation and extraction processes
* Spraying
* Digital-jet printing technologies

The overall dynamics of drop formation
--------------------------------------

* Mode 1: **Dripping** - at low efflux velocity - where liquid leaves the tube in periodic drops
* Mode 2: **Jetting** - at high efflux velocity - a coherent laminar jet emerges from the tube before breaking into discrete drops due to the **Rayleigh instability**

The dynamics of dripping
------------------------

* Stage 1: **Quasi-static growth of the drop** - ends when surface tension is unable to counter the weight of the drop
* Stage 2: **Necking and breakup:** 

- Caused by destabilising effect of surface tension and leads drop to develop a conical upper portion and almost spherical lower portion.

- The liquid thread connecting these two simultaneously increases in length and decreases in diameter to a **critical breaking point** at which the spherical portion is released to form the primary drop. 

- The remaining liquid column recoils rapidly due to surface tension and may then break up again at the top, leading to the formation of a **satellite drop**

**Q: If the industrial application involved any time dependency - the critical breaking point may be important**

**Q: Satellite drops may be important for accuracy? - or are these an extreme case?** 

Experimental studies
--------------------

Experimental study of Peregrine (1990) showed four stages:

* necking
* bifurcation
* recoil
* secondary necking

Shi (1994) showed increasing drop viscosity leads to formation of long thin threads - showing a microthread or a series of smaller necks can by spawned near breakup.

Zhang and Basaran (1995) showed effect of capillary radius and fluid properties on primary and satellite drops, for cases with and without surfactants - comparing water with 85% glycerine at low flowrate.

Modelling
---------

Macroscopic force balance
~~~~~~~~~~~~~~~~~~~~~~~~~

Early models attempted to predict drop size as a function of fluid properties, nozzle geometry and flowrate  - based on **simple macroscopic force balance** - where necking and breakup are preceded by pure static growth from the nozzle. However such models are inaccurate by >20%

One-dimensional modelling
~~~~~~~~~~~~~~~~~~~~~~~~~

Based on Rayleigh (1878) drop breakup was analysed by a range of 1D axi-symmetric models, this can provide information on **interface rupture during drop formation**. The method is to solve the axial mass and momentum conservation equations. These can predict the **time evolution of the drop shape** at low flowrates - showing how thin liquid threads above the primary drop can spawn a series of smaller necks with ever thinner diameter prior to breakup.

No 1D model can capture simultaneously:

1) Macroscopic features, e.g. the primary drop volume and time to breakup  
2) Microscopic features, e.g. thin liquid threads and satellite drops

Two-dimensional modelling
~~~~~~~~~~~~~~~~~~~~~~~~~

Boundary integral method / boundary element method of Schulkes (1994):

* Inviscid, irrotational or Stokes flow
* 2D axisymmetric

VOF method of Richards (1995): 

* liquid jet formation from a nozzle into a second liquid

Full axi-symmetric Navier-Stokes equations method of Wilkes (1999):

* Formation of a drop of Newtonian fluid at the finite Reynolds number encountered in ink-jet printing.
* Used FE method - allowing body fitted grid.
* Can simulate drop volume and thin filaments successfully.
* Chen (2002) showed that this can be used to simulate the **pinch off** of a low viscosity fluid where the surface overturns before pinch off.
 
VOF common because of the ability to simulate micro/macroscopic features.

Engineers don't have time or expertise to develop codes for complex free surface flows - so this study compares CFX and FLOW3D.

Method of comparision:

* Previous experiments
* Previous computational studies
* New experiments

In the end there is a compromise between:

* Computing time
* Accuracy

Experimental Method
===================

Drop formation visualised with high-speed photography with image analysis - quantitative information about:

* limiting length of the drop
* primary drop volume
* time to break-up

Apparatus
---------

* Constant flow rate to a capillary tube supplied by **syringe pump** 
* Capillary tube is on a vibration isolation table
* Ambient temperature is monitored using a thermometer
* Three sizes of capillary tube 0.4mm, 0.8mm, 1.375mm (outer radii)
* Square edge to tube ensures interfacial surface is at the outer edge.
* Length to diameter ratio >10 for fully developed flow at outlet
* Ratio of inner to outer radii > 0.3 :math:`\Rightarrow` effect of capillary wall thickness is negligible
* Pump produced flowrate in range 0.01 to 25 ml/min

**Really need dimensionless numbers here - e.g. what is the Reynolds number & Stokes number range?**

* Images captured with a 4500frame/sec camera (256x256 pixels) and 40500 frames/sec (64x64 pixels)
* Lighting with a 250W halogen light reflected off white card
* Use of a video monitor to view the images

Materials and Property Measurements
-----------------------------------

* Glycerine-water was used because the viscosity can be made to vary by three orders of magnitude as the concentration increases from 0 to 100, whilst the surface tension varies only by 10% and the density by 20%.

* Density of mixtures measured with a density bottle
* Viscosities measured using a rheometer
* Surface tension measured with a ring probe mounted on a torsion balance.

Procedure
---------

* Air is cleared and steady flow established.
* Drop formation at constant flow rate is extremely repeatable - volume at break off varies by only 5%

CFD Analysis of Drop Formation
==============================

Assumptions:

* Incompressible
* Newtonian
* Axisymmetric

Dimensional analysis:

* Used dimensionless velocity, lengths, stresses and time

Physics:

* Momentum Equation contains Reynolds number and Stokes Number
* Surface tension introduced via free surface boundary condition (i.e. viscous and surface tension stresses should be in balance)
* Applied a fully developed velocity profile at inlet
* No mass transfer across free surface

Mathematical description of free surface:

* Eulerian - fixed mesh - fluid moves through mesh
* Lagrangian - moving mesh - fluid convected with mesh

Some problem with Lagrangian approach - unable to predict drop dynamics past the point of necking and subsequent drop pinch off

Eulerian approach with VOF method tracks interface using a fluid marker function convected with flow. Fluid marked as 0, air as 1. Between 0 and 1 is the free surface. This requires a finer mesh than Lagrangian approach.

CFX - uses two phase approach - flow equations are solved in liquid and air regions

FLOW3D - uses single phase approach - free surface is a discontinuity - flow equations are only solved for liquid. Hence FLOW3D should be more efficient than CFX.

Results
=======

Experimental
------------

* Plotted time to breakup vs flowrate on log scale for different concentrations of glycerine.
* Agreement with previous results is excellent
* Viscous effects cause:
  
   - Elongation of liquid cone above primary drop
   - Lengthening of liquid thread

* For a fixed capillary size - decrease in time between breakup is 10% when viscosity is reduced from 50% glycerine to 20% glycerine

* Effect of increasing flow rate is a reduction in time between breakup - which decays exponentially, (except for smallest diameter tubes, which approaches jetting at higher flowrates)

* For small tube diameters, as flow is increased the inertia of the fluid forces more fluid into the droplet, increasing it's size. At larger flowrates, this inertia pushes the droplet further from the tube and ultimately a jet is formed.

* The droplet formation is now dominated by the Rayleigh instability occurring along a column of fluid, with a diameter of the order of the inner tube diameter, as opposed to the diameter of the order of the droplet, hence drop volume is reduced.

* **Important industrial driver is need to increase resolution by reducing drop size by reducing capillary size. However, we must avoid blockage** Ways of reducing drop diameter:

1) There is a diminished return on the reduction of drop volume by reducing tube diameter 
2) Control waveform used to drive a piezo element
3) Use a non-circular cross section for the tube

Numerical Simulation
--------------------

Initial drop formation
~~~~~~~~~~~~~~~~~~~~~~

* Big issue - smearing of the free surface at the interface
* CFX (multi-phase solution) shows slightly more diffuse results than FLOW3D (single phase solution)
* Marker cell function (or volume fraction) of 0.5 was best representative of the free surface

Drop development
~~~~~~~~~~~~~~~~

* CFX shows better agreement in terms of shape with previous simulations
* FLOW3D shows better agreement in terms of time taken to reach a drop shape with previous simulations
* FLOW3D is much more efficient than CFD (2000 versus 433000 CPU seconds)

Comparison with visualised drops - 20% glycerine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Time between drop breakup shows good comparison
* Agreement deteriorates as drop break up approaches

Comparison with visualised drops - 50% and 20% glycerine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Higher flow rates and higher viscosity liquid mean than the thread is not present - where surface tension isn't dominant
* When surface tension is dominant, the comparison is poorer

**Failure may be due to CSF algorithm used by FLOW3D, which represents surface tension via a body force, which is known to have limitations for surface tension-dominated flows.**

Finite element schemes are able to capture this long neck feature

Comparison of drop length
~~~~~~~~~~~~~~~~~~~~~~~~~

* Drop length versus % glycerine shows a good comparison at low viscosity
* At high viscosity the error is larger
* Unable to obtain converged solutions for very low viscosity - i.e. water

Comparison of drop volume
~~~~~~~~~~~~~~~~~~~~~~~~~

* Generally good comparison, with some differences

Comparison of time to breakup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Time to breakup is in excellent agreement between simulation and experiment

Prediction of satellite drop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* This is predicted, however, requires a much finer grid - probably requires local grid refinement

Conclusions
===========

Drop formation dynamics is important because:

* It has a large array of industrial applications
* The nature of the free surface flow is complex

Justification of approach:

* There are accurate methods for simulating drop formation, but they are too complex for engineering applications
* Hence this study compared CFX and FLOW3D

Physics:

* Drop length is sensitive to viscosity at higher viscosities
* Drop volumes are less sensitive to viscosity
* Drop volume has a complex relation with flow rate
* Drop volume can be reduced by reducing tube diameter, but at a diminished return

Prediction using CFD:

* At higher flowrates and viscosities CFD predicts well (**due to greater influence of viscous effects over surface tension?**)
* It can predict time to break-up, primary drop volume
* Not good at predicting interfacial shapes, due to smearing at free surface (**due to volume fraction averaging?**)
* It can't predict long thin threads between liquid and cones at low flows and satellite drops on breakup
* More generally CFD packages perform poorly where there are low flowrates and high surface tension forces

Efficiency:

* FLOW3D is much more efficient than CFX, because of the use of a single phase with a discontinuity (FLOW3D) versus two phases (CFX)


