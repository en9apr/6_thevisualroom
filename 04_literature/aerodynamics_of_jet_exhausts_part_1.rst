===================================
Aerodynamics of Jet Exhausts Part 1
===================================

Goulos, I., Stankowski, T., Otter, J., MacManus, D., Grech, N. and Sheaf, C. (2016) ''Aerodynamic Design of Separate-Jet Exhausts for Future Civil Aero-engines - Part 1: Parametric Geometry Definition and Computational Fluid Dynamics Approach`` Journal of Engineering for Gas Turbines and Power, Vol 138.

.. contents::
   :local:

.. highlight:: latex

Abstract
========

- Output is: **integrated approach** targeting **aerodynamic design** of **separate-jet exhaust systems** for future gas-turbine aero-engines.
- Framework is a series of fundamental theories applicable to:
    * engine performance simulation
    * parametric geometry definition
    * viscous/compressible flow solution
    * **design space exploration (DSE)**
- Method:
    * Mathematical method developed based on:
        + **class-shape transformation (CST)** functions for geometric design of axi-symmetric engines
        + Standard set of nozzle design parameters
    * Design carried out using:     
        + Flow capacities established from **0D cycle analysis**
    * Coupled to:
        + ICEM for automatic mesh generation using **block structured approach**    
        + Fluent for **RANS solution**
    * Validation against:
        + Experimental data on a small-scale turbine powered simulator (TPS)
    * Coupled tool to:
        + **DSE Latin-hypercube sampling**
    * Applied to two civil engines:
        + Current
        + Future
- Results:
    * Relation between **exhaust systems thrust** and **discharge coefficient** has been quantified
    * **Dominant design variables** that affect aerodynamic performance of the exhausts have been determined
    * Comparative evaluation of the **optimised exhaust design** of each engine
- Conclusions:
    * Enables aerodynamic design of exhausts using only a **few design variables**
    * Enables **quantification and correlation** of aerodynamic behaviour of each engine architecture
    * Is an enabling technology to **identify fundamental aerodynamic mechanisms** for exhaust system performance
    
Introduction
============

Background
----------
- What is the future trend in civil turbofans?
    * The **motor** of civil turbofan engines will have greater **thermal efficiency**:
        + Increased TET
        + Increased OPR
    * Maybe leads to intercooled and intercooled-recuperated cycles
    * Future turbofan engines will have lower **specific thrust** and improved **propulsive efficiency**:
        + Higher BPR (15+, it is currently ~11)
        + Lower FPR
- Why is the **exhaust** important?
    * Higher BPR means higher gross to net propulsive force ratio
    * High BPR designs are therefore more sensitive to gross propulsive thrust
    * Gross propulsive thrust is linearly dependent on the aerodynamics of the exhaust
- Why is the **bypass** duct important?
    * High BPR means higher mass flow through bypass
- **Post-exit components** are also important

Performance Prediction of Engine Exhaust Systems
------------------------------------------------

- Engine housing is not designed by engine manufacturer, so **thrust-drag bookkeeping** (TDB) is needed to mitigate losses.
- Exhaust system can cause 1.5 to 2\% loss in gross propulsive thrust
- In TDB :math:`C_V` (velocity coefficient) and :math:`C_D` (drag coefficient) are used for measuring performance
- CFD used for aerodynamics analysis of exhaust nozzles
- What are the flow features?
    * Boundary and shear layer interaction
    * Expansion waves
    * Shock waves
- What is the accuracy of CFD?
    * less than 1% for :math:`C_D` and :math:`C_V`, largely due to uncertainty in exprimental data

Scope of Present Work
---------------------

- What is unique about the current work?
    * Methodological approach for:
        + Parametric geometry definition
        + Aerodynamic analysis
        + Examination of separate jet exhaust systems
    * Impact of high BPR and lower FPR on exhaust system design
    * **Not considering installation geometry then?**
    
- What are the objectives of the current work?

    * Derive analytical formula for **parametric geometry definition** of separate jet exhausts
    * **CFD model** of bypass duct, nozzle and post exit conditions
    * **Framework for exploring design space** for aerodynamic performance
    * Explore design space for **future and current engines**
    
- How is the parametric geometry defined?

    * CST functions (class function shape function transformation)
    * Axi-symmetric
    * Separate jet exhausts
    * Extends Qin's aerofoil approach to exhausts and nozzles
    * Parameterisation based on required flow capacities
    * Coupled to ICEM and Fluent

- How is the CFD model defined?

    * CFD validated against small scale turbine power simulator (TPS)
    * **What is the definition of the CFD model? (section below)**
    * **BCs, discretisation scheme, solver, turbulence model?**

- How is the design space defined?

    * Coupled to framework
    * Explores future and current turbofan
    * **How is DSE done (second paper)?**
    
Numerical Approach
==================

Methodological Overview
-----------------------

- What is GEMINI?

    * **Geometric Engine Modeler Including Nozzle Installation**
    * Designs separate jet exhaust systems based on key **engine hard points**
    * Applicable to:
        + Engine performance simulation
        + Exhaust nozzle geometry
        + Parameterisation
        + Viscous compressible flow solution

- How is the 0D engine performance model defined?

    * Inputs: thermodynamic and geometric design parameters
    * Analyse engine cycle at design point and off design
    * Uses Cranfield's Turbomatch
    * Outputs: **size** of bypass and core, **average flow properties** at inlet and exit of bypass and core

- How is the GEMINI, ICEM, Fluent and Post processing done?

    * Inputs: **flow capacities** and **size** of bypass and core
    * Inverse design approach in Gemini produces 2D axi-symmetric geometry
    * Transfers to ICEM
    * Transfer to Fluent
    * Transfer to Post processor
    * Outputs: :math:`C_D^{bypass}` and :math:`C_D^{core}` and :math:`C_V^{overall}`

Engine Performance Simulation (Turbomatch)
------------------------------------------
    
- How is the 0D engine performance model done?

    * Turbomatch
    * 0D aerothermal analysis
    * Solves for mass and energy balance between engine components
    * Assumes engine is operating at steady state
    
Parametric Geometry Definition of Exhaust Nozzles
-------------------------------------------------
    
- How is the parametric geometry defined?

    * Kulfans CST functions
    * Qins CST (class shape transformations) extended from aerofoils to exhausts
    * nth order Bernstein polynomial - **uses a summation of polynomials** to describe the surface **with an offset for position**
    * The geometry is split into the **upstream duct** and **exhaust nozzle**
    * Geometric parameters are specified to achieve design parameters using **control points** (where geometric information is avaliable)
    * :math:`(n-1) \times (n-1)` system of linear equations created
    * BCs are established from control points
    * **How is the geometric BCs satisfied to be unique? (e.g. is the gradient specified as well?)**
 
CFD Domain and BCs
------------------

- 2D axi-symmetric
- Why is the engine intake included?
    * Domain includes engine intake to account for effect of mass flow capture ratio on the nacelle pressure distribution
    * This is required to capture the static pressure aft of the nacelle afterbody and the effect of freestream supression on the aerodynamics
- Freestream:
    * Pressure far-field
    * static pressure, static temperature, Mach number
    * Position of freestream: 150 maximum nacelle diameters **Is this really big enough?** (despite sensitivity analysis, maybe ok if inviscid)
- Fan face:
    * Pressure outlet
- Bypass:
    * Pressure inlet
- Core:
    * Pressure inlet
- Vent:
    * Prescribed mass flow
- How is the non-uniformity of flow accounted for?
    * **Streamline curvature method** applied to fan rotor and fan outlet guide vanes

Automatic mesh generation
-------------------------

- Block-structured mesh automatically generated using ICEM
- y+ is unity
- 50 nodes normal to aeroline surface
- Expansion ratio 1.2
- Mesh topology based on MSc thesis?
- **Why not use more efficient hybrid mesh generation?** 
- **Why not use better scripting language than ICEM e.g. Pointwise?**
- **Why not use better quality expansion using hyperbolic PDE in boundary layer using Pointwise?**

Definition of CFD Approach
--------------------------

- ANSYS Fluent
- RANS using :math:`k-\omega` SST turbulence model
- Green-Gauss for gradients
- 2nd order upwind scheme for flow variables, turbulent kinetic energy and dissipation rate 
- Thermal conductivity via kinetic theory
- Eighth order polynomial for specific heat capacity (:math:`C_P`)
- Sutherlands law for dynamic viscosity
- **Why not MUSCL scheme?**
- **Why not Riemann solver instead of slow SIMPLE algorithm?**
- **Acoustics cannot be included using steady state CFD model**
- **Solution won't be solver independent**

Exhaust System Performance Accounting
-------------------------------------

- Discharge coefficient:

.. math::

    C_D = {{\dot{m}_{actual} }\over {\left( {\dot{m} \over A }\right)_{ideal} \ A_{throat}}}

- The throat area is taken to be equal to the exit area **Is this valid? Is there a vena contracta?**
- **It could be like a Venturi meter, where the contraction coefficient is unity, such that** :math:`C_D` **equals** :math:`C_V` **a ratio of velocities for single phase flow**
- :math:`C_D` is defined for the core and the bypass separately

- Gross propulsive force:

.. math::

    F_G = F_G^{bypass} + F_G^{core}  + F_G^{zone3} - \text{integral of (static pressure term in axial direction aft of max nacelle diameter - viscous shear stress)} 

- Overall velocity coefficient (divide a force by a mass flow rate and you get the actual velocity on top):

.. math::
    
    C_V^{overall} = {F_G \over { \left(\dot{m}_{actual}^{bypass} V_{ideal}^{bypass} + \dot{m}_{actual}^{core} V_{ideal}^{core} + \dot{m}_{actual}^{zone3} V_{ideal}^{zone3} \right) }}
    
Results and Discussion
======================

Grid Sensitivity Analysis
-------------------------

- Numerical predictions at DP mid cruise conditions
- 5 meshes using uniform refinement
- Around 100,000 cells for coarse mesh, 1 million cells for fine mesh
- **Non-monotonic behaviour could be caused by turbulence model**
- **Non-montone behaviour due to limiter in 2nd order scheme?**
- **Investigate the effect of higher order schemes on monotonicity?**
- **May be able to use coarser mesh with 3rd order scheme?**
- **Big Problem: No AMR - may be able to use even very coarse grid with AMR and high order scheme**

Validation of Employed CFD Approach
-----------------------------------

- Pylon blockage in experiment, so CFD must be corrected
- No correction for 3D nature of flow, CFD is 2D axi-symmetric
- Used different FPRs and measured normalised mass flow and gross propulsive thrust
- Difference is around 5\% due to 3D nature of flow and possibly uncertainty about pylon
- Isentropic Mach number is around 10\% different in bypass and 6\% in core
- **Possibly because of lack of resolution around shock waves?**

Design Space Exploration
------------------------

- Design of Experiment approach is **Latin Hypercube** to mitigate the cost of CFD simulations
- After a representative database is collected, the beaviour is investigated statistically
- Design variables are correlated with the performance metrics using **Pearson's product moment of correlation**

Case Study Description
~~~~~~~~~~~~~~~~~~~~~~

- Two engines, Current (E2) and Future (E1) with BPR of 11 and 16 respectively
- Each cycle has been optimised wrt FPR to maximise **specific thrust** and minimise **specific fuel consumption**
- **How was it optimised?**
- DP mid cruise conditions for both engine models
- Bypass is choked, core is unchoked

Design Space Definition
~~~~~~~~~~~~~~~~~~~~~~~

- 11 and 12 parameters for E1 and E2 engines have specified ranges, in agreement with **design guidelines** and **manufacturing constraints**

Preliminary Statistical Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Each design space was discretised using the Latin Hypercube method
- 360 exhaust geometries were used per engine
- Correlation between imposed design variables and performance metrics was investigated
- Question: Which are the dominant variables?
    
- Large percentage variation in **core discharge coefficient** and **zone 3 pressure ratio**, due to **strong influence of core cowl design on core nozzle exit static pressure**    
- E2 has an additional parameter, giving it **more degrees of freedom** than E1, so the variation in the values is greater
- Definition of velocity coefficient renders it relatively independent of discharge coefficient to first order, leading to smaller standard deviation for the velocity coefficient.
- **Why did E2 have more degrees of freedom?**

Assessment of Apparent Design Space Linearity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Plotted charts and determined Pearson correlation coefficient for:

    * :math:`C_V^{overall}` versus :math:`C_D^{bypass}`
    * :math:`F_N` versus :math:`C_D^{bypass}`
    * :math:`F_N` versus :math:`C_V^{overall}`
    
- Exchange rates between :math:`F_N` and :math:`C_V^{overall}` can be almost double for future engines compared to current engines
- Hinton Diagrams for all performance metrics versus all design variables, coefficients are dependent only on three main design variables
- Increasing nozzle :math:`C_P` to exit length ratio **moves low pressure turbine hump upstream** and mitigates **strong shock**
- This improves discharge coefficient by 0.4 \% and velocity coefficient by 0.06\% and increases :math:`F_G` by 0.45\%

- **Why are the improvements so small?** But I suppose nearly 0.5\% is large for discharge coefficient?

Conclusions
===========

- Integrated approach for aerodynamic design of separate jet exhaust systems
- Applicable to:
    * Engine performance simulation
    * Parametric geometry definition
    * Viscous compressible flow
- Analytical approach for parametric geometry using CST functions
- Validated against experimental data
- Formulation for design space evaluation
- Used future and current aero engines
- Sensitivity to parametric changes has been identified
- Hinton diagrams are effective in representing behaviour and to identify guidelines for design
- Can be used to identify fundamental aerodynamic mechanisms




