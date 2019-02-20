======================================================================
Wen and Mobbs (2014) "Numerical Simulations of Laminar Air-Water Flow"
======================================================================

Wen, X., Mobbs, S., "Numerical Simulations of Laminar Air-Water Flow of a Non-linear Progressive Wave at Low Wind Speed", Boundary Layer Meteorology, 150 (2014) pp 381 - 398

.. contents::
   :local:

Where is this study placed within the wave modelling community?
===============================================================

* This paper is looking at the the issue of the generation of waves by wind. Unlike the Miles (1957) formulation and the subsequent development, it allows non-linear effects. It is known that non-linear effects are important when the phase speed and wind speed match (at the critical height), which is one of the cases reported here. 
* It is shown that the maximum wind speed moves closer to the surface as the wind speed increases, similar to the non-separated sheltering mechanism from the work of Jeffreys (1924). This has been shown by Belcher and Hunt (1998), in that for c/U < 15 (slow moving waves) and fast waves the mechanism for momentum transfer is non-separated sheltering.
* Flow regimes are defined based on wind speed :math:`U = 0` (no wind), :math:`U = u_m`,  :math:`U = c` to define in terms of the velocity field how momentum is transferred - this has not been defined before.
* Although the model includes a beach slope, it is not concerned with dissipation in shallow water, because an artificial viscosity is applied to damp out wave energy in the beach region (to allow the derivative of velocity to be zero at the outlet) and it also considers a deep water case (h/L = 0.6).
* Tsimring (1983) showed that it is valid to analyse components separately because the effect of other components is proportional to the air-sea density ratio times the square of the wave spectrum.

Summary of Findings
===================

* In the water, the potential flow versus viscous flow showed good agreement
* In the air, potential flow is different to viscous flow because of the effects of viscosity
* 3 wind speeds: :math:`U = 0, U = u_m` and :math:`U = c` are needed to compute the three different flow patterns

Literature Review
=================

* Studying air-water flows is important for the coupling between the atmospheric boundary layer and the oceanic boundary layer

What experimental measurements have been made for air and water oceanic surface waves?
--------------------------------------------------------------------------------------

* Suspended particle technique - velocity in air
* Growth of water waves?
* Turbulent velocity distribution in water under air-water interface
* Pressure in air at two fixed points above water level
* PIV - velocity distribution beneath wind driven air-water interface and analysed tangential stress
* Wave growth and attenuation in **following and opposing** airflow
* PIV - velocity in air and shear stress
* Wind induced growth of slow water waves
* Coupled boundary layer air-sea transfer (field campaigns)

What are the difficulties?
--------------------------

* Measurements can only be made well above the peak of the wave - hence lose surface detail
* Velocity distribution, pressure and shear stress distribution of water and air **near surface** is lacking

What theoretical models have been produced?
-------------------------------------------

* Wave generation theory (air-sea flow)
* Perturbation solutions of airflow equations - airflow over a wavy surface
* Ocean wave models - energy balance equations uncertainty due to empirical calculation of energy input, wave-wave interaction and energy dissipation - because mechanism of wave-air interaction not understood

What CFD models have been produced for air flows?
-------------------------------------------------

* Boundary Integral Method for potential flow
* Numerical simulations for viscous flow over a stationary wavy boundary - airflow over a travelling wave isn't the same as airflow over a stationary wave.
* Turbulent flow over a small wave - 2 equation model
* Effect of wave steepness on flow pattern in air above Stokes wave employing 2 equation model
* Low Re model for air flow over a sinusoidal wave
* Low and high order turbulence models over a Stokes wave
* DNS for flow over a sinusoidal wave
* LES for airflow following and opposing a sinusoidal wave

How do CFD models help?
-----------------------

* Help to determine air-water flwo in the vicinity of the interface
* Shape and movement of wave surface described by the analytical solution to Laplace's equation for potential flow

How are CFD models solved?
--------------------------

* Separate solution for air-water
* Simultaneous solution for air-water

What is the result for low wind speed?
--------------------------------------

* Upward momentum transfer

Aim
---

* Investigate effect of wind speed on air-water flow pattern for a non-linear progressive wave

Method
======

* As for Wen (2012)
* Damping zone introduced for artificial viscosity at beach to shorten computational domain

Results
=======

**Linear vs Non-linear:**

* Non-linear solution shows sharper peaks and flatter troughs
* Good agreement between Fenton's fifth order solution and the non-linear viscous wave
* Linear solution at inlet quickly evolves into non-linear solution as wave propagates - but use fifth wave from inlet

**For potential flow:**

* Discontinuity in the velocity field - due to the absence of viscosity
* Excellent agreement between potential flow and viscous flow <1%
* Differences occur at surface of water

**For viscous flow:**

* All air moves with positive velocity
* As wind speed increases, vortices roll forward
* Orbiting water has a strong effect on the air

**Zero wind speed:**

* Vertical and horizontal velocity almost equal when wind is zero
* Strongest upward momentum transfer at zero wind speed

**Wind speed less than maximum orbital velocity:**

* Two recirculations - one above crest, one above trough

**Wind speed greater than maximum orbital velocity:**

* Recirculation above crest disappears
* Air above crest less than maximum orbital velocity

**Wind speed equals phase speed:**

* Air above crest equals maximum orbital velocity

**Wind speed greater than phase speed:**

* Air above crest has greater velocity than max orbital velocity
* As the wind speed increases the maximum wind-speed moves closer to the surface of the water - like the acceleration of air over the tops of hills






