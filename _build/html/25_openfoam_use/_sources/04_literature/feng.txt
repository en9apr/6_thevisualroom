========================================================================
Feng at al. (2011) "Evaluation of Cloud Convection and Tracer Transport"
========================================================================

Feng W., Chipperfield M.P., Dhomse S., Monge-Sanz B.M., Yang X., Zhang K., Ramonet M. "Evaluation of Cloud Convection and Tracer Transport in a Three-dimensional Chemical Transport Model" Atmospheric Chemistry and Physics, vol 11 pp 5783-5803, 2011

.. contents::
   :local:

Where is this study placed within the wave modelling community?
===============================================================

* This study is really outside the wave modelling community, looking at meteorological processes
* However, it might take input from wave models such as the pressure drag (or form drag) and shear stress/viscous drag (or friction drag)

Aim
===

Diagnose the following:

* the updraft mass flux
* convective precipitation
* cloud top height

Compare CTM with three analyses:

* ERA-40
* ECMWF Operational
* ECMWF Interim

Method
======

Cumulus parametrisations:

* Convective adjustment
* Mass-flux schemes

Large number of processes involved:

* Chemistry
* Photolysis
* Aerosol
* Large-scale advection
* Convection
* Dry/wet deposition
* Planetary boundary layer mixing
* Emissions

CTMs use meteorological data. But the key uncertainty in tropospheric CTMs is:

* Accuracy of sub-grid scale transport by convection

Method is to:

* Use archived mass fluxes
* Assess impact of resolution
* Try different external forcing meteorology and surface data
* Try different parametrisations
* Compare TOMCAT with ECMWF reanalyses and by using Radon as a model tracer

**Note on wind modelling (as this may related to wind-wave interaction):**

* TOMCAT reads winds as spectral coefficients of vorticity and divergence. These are then averaged onto whatever model grid is being used as part of the spectral transform. If the forcing winds are higher resolution than the model grid then information from the higher wavenumbers is not used - the spectral coefficients are truncated

Datasets
========

* Convective mass flux
* Cloud top height
* Convective precipitation
* Radon measurements

Conclusions
===========

**Convective mass flux**

* Use of archived mass fluxes improves CTM
* TOMCAT CTM underestimates convective mass flux compared to archive
* Most severe disagreement is in vertical extent of convection. Archive shows strong convective transport up to 100hPa, but model extends to only 200hPa (1hectoPa = 1millibar)
* Resolution of CTM didn't make much difference to the convection

**Cloud top height**

* Cloud top height - not as greater difference as convective mass flux

**Convective precipitation**

* Model of convective precipitation generally captures latitudinal precipitation

**Radon**

* Resolution did make a difference to the Radon results
* Agreement up to 10km in middle latitudes


