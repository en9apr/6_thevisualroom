===================================
Aerodynamics of Jet Exhausts Part 2
===================================

.. contents::
   :local:

.. highlight:: latex

Introduction
============

Background
----------

- Civil aviation traffic will increase in the future, so reduce:
    * direct operating costs
    * fuel burn
    * emissions
    * noise
- Need to improve design technology in:
    * Motor
    * Propulsor
- Improve motor thermal efficiency by:
    * Increase TET
    * Increase OPR
- Future architectures:
    * **Higher BPR** (currently 11, will be 15+)
    * Lower FPR
- This will:
    * **Lower specific thrust**
    * **Improve propulsive efficiency**
- Why is the **exhaust** important?
    * Increasing **BPR** will increase **gross to net propulsive thrust**
    * Designs are therefore more **sensitive to variations in gross propulsive thrust**
    * Gross propulsive thrust is **linearly dependent** on the **aerodynamic performance of the exhaust**
- Which components will be analysed?
    * **Bypass duct**
    * **Nozzle**
    * **Post exit components**

Exhaust system performance accounting
-------------------------------------

- What is separate jet exhausts?
    * **Core cowl** separates core flow from bypass flow
    * **Protruding core plug** limits the length of core cowl
- What is the problem with separate jet exahusts?
    * Can be substantial sources of thrust loss (gross thrust reduced by 2\%)
- How is performance measured?
    * Discarge coefficient
    * Velocity coefficient
    
Design Optimisation of Engine Exhaust Systems
---------------------------------------------

- CFD is a reliable performance prediction tool
- It is also efficient and used in design optimisation (**Efficiency is depending on mesh size, geometry complexity, schemes, solvers and models used**)

Heath (2015)
~~~~~~~~~~~~

- Axi-symmetric, dual stream plug nozzle
- Parametric geometry via free-form deformation and 3rd order b-splines
- RANS solver
- **Steady state**
- Unstructured grid
- Adjoint, grid deformation, grid adaption to obtain gradients
- Optimisation using sequential quadratic programming (SQP)
- Minimise integral of **near-field pressure disturbances** relative to freestream flow
- Gross thrust gain of 0.2\% relative to baseline

Clemen (2012)
~~~~~~~~~~~~~

- **Why not use HYDRA for?**
    * Linearised unsteady solver
    * Non-linear solver
    * Steady adjoint solver
    * Harmonic adjoint solver
- Integrated framework for high BPR turbofan with core mounted gearbox **Like Ultrafan?**
- 2nd order splines for parametric geometry
- **3D RANS solver HYDRA**
- **Steady state**
- Hybrid optimisation comprising initial design of experiment coupled with RSM and global optimiser
- **RSM** (Response Surface Modelling) based on design of experiment results using interpolations based on radial basis functions
- **Genetic algorithm** to minimise total pressure loss within bypass duct
- 0.1 \% reduction in pressure loss

Haderlie and Crossley (2010)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Axi-symmetric supersonic inlet
- Modified splitter geometry that separates core and bypass flow
- Parametric geometry based on Kulfan's CST method
- RANS flow field 
- Multiblock structured mesh 
- Optimisation based on design of experiment - Latin hypercube
- Surrogate model using Kriging interpolation
- Optimisation from genetic algorithm and local gradient based sequential quadratic programming
- Optimisations used total pressure recovery and peak radial distortion intensity at the inlet's aerodynamic interface plane
- Improved splitter design that satisfied imposed geometric constraints
- **Current paper is based mainly on this one**

Qiu (2014)
~~~~~~~~~~

- Unsteady, continous adjoint-based acoustic propagation method
- Optimise the design of a low bypass duct for a civil turbofan
- Hick-Henne shape functions for parametric model of bypass and nozzle
- Optimisation based on local gradient based algorithm driven by Jacobian from adjoint method
- Minimise tonal noise 
- Reduced overall SPL in far-field by 2.78dB

Scope of Present Work
---------------------

- Aerodynamics of the exhaust is important for future high BPR engines
- What is unique about the current work?
    * Previous authors have looked at optimising exhaust nozzles
    * A holistic approach for **separate jet exhausts** including bypass, core duct and post exit components has not been reported
    * Impact of **high BPR engines** and lower FPR on exahust system design and optimisation has not been reported
- What is the approach?
    * Cycle analysis
    * Geometry parameterisation
    * Mesh generation
    * RANS flow solution
- What is new?
    * Expand optimisation strategy using DOE (Design of Experiment), RSM (Response Surface Modelling) and GA (Genetic Algorithm)
- What is being optimised?
    * Current and future engine architectures
    * Large turbofans
    * Optimise the exhaust designs

Numerical Approach
==================

Aerodynamic Design of Separate-Jet Exhausts
-------------------------------------------

- What is GEMINI?
    * Tool developed is GEMINI
    * Designs complete exhaust system for designated engine cycle using key engine hardpoints
    * Applicable to:
        + Engine performance simulation
        + Exhaust duct and nozzle aeroline parameterisation
        + Viscous compressible flow
- What is the process?
    * Designate a set of thermodynamic cycles and geometric design parameters
    * Analyse engine at design point and off design (0D conditions) - Turbomatch, output bypass and core sizes and flow capacities, at steady state conditions
    * Inverse design to create 2D axi-symmetric model
    * Automatic generation of grid
    * Convege CFD solution
    * Determine discharge and velocity coefficients

Exhaust System Parametric Geometry Definition
---------------------------------------------

- How is the parametric geometry defined?
    * Kulfan's CST functions
    * Qin's CST variations
    * Bypass, core, duct exhaust are reduced to a set of analytical expressions
    * The expressions are functions of a standard set of design parameters

- How is the nozzle designed?
    * Geometric throat area is known
    * An effective convegent-divergent ratio is defined
    * Application of the rolling ball area estimation method to nozzle exit plane and upstream CP results in a series of control points that satisfy the prescribed design parameters
    
- How is the upstream duct defined?
    * Direct control of a series of control points
    
- Why is the engine intake considered?
    * To capture the effect of inlet mass flow capture ratio
    * To then account for the effect of the static pressure distributionon the nacelle
    * To then account for the effect of freestream supression on the aerodynamic performance of the exhaust system

- How is the geometry defined?
    * Upstream duct via specifying position, slope and curvature within a series of control points
    * Core cowl and plug are modelled as straight lines
    * Includes a third nozzle

DSE and Optimisation
--------------------

- What is done in this paper?
    * Extend GEMINI
    * Implement DSE and optimisation environment
    * Non-linear nature must be dealt with
    * Must mitigate the cost of numerous CFD applications
    
- How is the process of DSE done?
    * Deployment of **DOE** method to explore the available design space
    * Construct **RSMs** from DOE results

- What kind of DOE is used?
    * **Latin Hypercube**

- What is a RSM?
    * Hypersurface describing the mathematical relationship between a set of imposed design inputs and outputs
    * The use of RSMs will avoid a **prohibatively large number of CFD simulations**
    * Interpolation using Gaussian process regression, **Kriging interpolation**
    * Performance metrics are **discharge and velocity coefficient**
    * **Leave-one-out cross validation** used to check predictive accuracy of RSMs
    
- How is the optimisation done?
    * Global method to avoid being trapped in locally optimal solution - GA (Genetic Algorithm)
    
Results and Discussion
======================

Definition of Baseline Engines
------------------------------

- How are the baseline engines defined?
    * Optimise low pressure exhaust system design and core afterbody aerolines for current and future aero-engines. 
    * BPR current = 11
    * BPR future = 15
    * OPR, TET, component efficiencies selected according to technology guidelines
    * Each cycle optimised wrt FPR to maximise specific thrust
    * 2D axi-symmetric
    * Geometry from public domain
    * Predictions at mid-cruise
    * Bypass is choked, core is unchoked
    
Parametric Design Space Definition
----------------------------------

- How are the parameters designed?
    * 11 to 12 variables for future and current engines
    * **Outer line angle** is kept constant for future engine
    
Design Space Exploration
------------------------

- How is the design space explored?

    * Design space discretised using Latin Hypercube
    * 360 exhaust geometries 
    * Correlation between design variables and performance metrics was investigated
    * Hinton Diagram using Pearson's product-moment correlation
    * Shows only a few parameters influence the performance metrics

- Within the range of assumptions, the aerodynamic **performance of the exhaust is decoupled from the intake and nacelle forebody**
- Changes applied to the exhaust do not influence the intake or nacelle

Response Surface Modelling
--------------------------

- How are the RSMs constructed?
    * Using DOE data
    * Interpolation using Gaussian processes regression (Kriging interpolation)
    * Quadratic regression function and squared exponential autocorrelation function
- How are the RSMs checked?
    * Leave one out cross validation
    * Employs all the avaliable data apart from one, which is the one to prediction
    * Prediction is compared with original raw data for accuracy
    * Surrogate model predictions are correlated against raw data using Pearson's product moment of correlation
    * Also assesses averge model error and standard deviation for each performance metric
    
- Result shows that CFD raw data and predicted data has very high correlation.
- **Could be improved using a larger amount of data**
- Low percentage error
- Standard deviation is of similar order to error, so data is scattered - shows non-linearity of the system

Exhaust System Design Optimisation
----------------------------------

- How is the optimisation performed?
    * Genetic algorithm
    * Advantage of using RSMs is that they are more efficient than CFD models
    
- What is the process for the GA?
    * For current and future engines
    * Optimise in terms of overall velocity coefficient
    * Population size is 10 times number of design variables
    * 40 generations
    * Convergence criterion of :math:`10^{-12}`
    
- What are the results of the optimisation?
    * Good solution achieved within 500 evaluations
    * Still contains small number of unfit individuals
    * Improvement wrt baseline values is large (2-4% in thrust)
    * CP to exit length ratio is increased (as before) mitigating strong shock
    * Also flow separations are mitigated
    
Conclusions
===========

- Design optimisation for separate jet exhausts for future civil aero engines
- Modules for:
    * Cycle analysis
    * Geometry parameterisation
    * Mesh generation
    * Viscous compressible flow solution
- Novel analytical geometry tool using CST functions
- 2D axi-symmetric RANS CFD model
- Extended formulation to include:
    * DOE
    * RSM
    * GA
- Used to optimise:
    * Current engine
    * Future engine
- Design optimisation can increase net propulsive force by 1.4\% or 3.4\% for future and current engines
- Can identify design guidelines and mitigate undesirable flow features








