============
OpenFOAM Use
============

This is the use of the OpenFOAM language.

.. contents::
   :local:

What is the process in CFD?
---------------------------

* Analyse problem 
    - literature survey
    - dimensionless groups
    - physical principles
* Model physics 
    - NSE, turbulence, other physics
    - Identify models and coefficients
* Setup geometry and mesh
* Determine numerical method
    - transient/steady state
    - differencing schemes
* Run calculation
* Analyse and post-process results

Directories
-----------
    

    
How to get to the installation directory?
=========================================

The following command gets you to /opt/openfoam4/

::

    $ foam    

    
How is the installation directory structured?
=============================================

Shortcut to ``/opt/openfoam4`` = ``$WM_PROJECT_DIR``

* /opt/openfoam4/ 
    - applications
        + solvers
            * basic
            * discreteMethods
            * financial
            * lagrangian
            * combustion
            * DNS
            * heatTransfer
            * multiphase
            * compressible
            * electromagetics
            * incompressible
            * stressAnalysis
        + test
        + utilities
    - doc
    - src
    - tutorials
    - bin
    - etc
    - wmake
    - platforms

    
How to get to the run directory?
================================

The following command gets you to /home/andrew/Dropbox/7_OpenFOAM/run/

::

    $ run  
    
How is the user directory structured?
=====================================

Each of 1_lid_driven_cavity etc represents a "case directory".

* /home/andrew/Dropbox/7_OpenFOAM/run/
    - 1_tutorials/
        + 1_lid_driven_cavity 
        + 2_pipe_flow
        + 3_flat_plate
* /home/andrew/Dropbox/7_OpenFOAM/applications/
* /home/andrew/Dropbox/7_OpenFOAM/lib/    
    
What is in the case directory?
==============================

A list of dictionaries:

* /constant/
    - turbulenceProperties
    - transportProperties
    - /polyMesh/
        + boundary
        + faces
        + neighbour
        + owner
        + points
* /system/
    - blockMeshDict
    - controlDict
    - fvSchemes
    - fvSolution
* /0/ (timesteps)
    - U
    - p
    - k
    - epsilon
    
What is the constant directory?
===============================

* turbulenceProperties - turbulence properties
* transportProperties - physical properties e.g. viscosity
* polyMesh - mesh details
    
What is in the system directory?
================================

* blockMeshDict - mesh scaling, basic meshing, boundary names
* controlDict - application, time start/end, timestep, write controls
* fvSchemes - differencing schemes
* fvSolution - solver, solver controls, pressure-velocity coupling controls

What is in the 0 (timestep) directory?
======================================

* Boundary conditions
* Initial conditions  

Dictionaries
------------

What are dictionaries?
======================

A dictionary is an entity that contains data entries that can be retrieved by the I/O by means of keywords. The keyword entries follow the general format (keyword-value pairs)

::

    <keyword>  <dataEntry1> … <dataEntryN>; 

* Free format ASCII text files
* Parsed by OpenFOAM
* Only values actually needed are read in
* Comment lines are ``//`` or ``/* .... */`` (C format comments)
* More flexible than alternatives
* Format:
    - banner
    - FoamFile - dictionary with general information
    - e.g. transportProperties - the dictionary with the actual information
    
What is a dimensionedScalar in the transportProperties dictionary?
==================================================================

A scalar with a:

* Name
* Physical dimensions
* Value

.. code-block:: c

    nu [0 2 -1 0 0 0 0] 0.01;
    
How are the dimensions given in the transportProperties dictionary?
===================================================================

.. code-block:: c

    [M L T theta mol I Cd]
    
* M = mass (kg)
* L = length (m)
* T = time (s)
* theta = temperature (K)
* mol = quantity (mol)
* I = current (A)
* Cd = luminous intensity (candela)

OpenFOAM checks dimensions of equations


What is the turbulenceProperties dictionary?
============================================

Read by pisoFoam and simpleFoam to switch turbulence model, e.g.

.. code-block:: c

    simulationType RASModel;
    
    RAS
    {
        RASModel kEpsilon;      // chosen model
        turbulence on;          // turbulence on
        printCoeffs on;         // coefficients printed to terminal
    }

Coefficients are stored in sub-dictionaries with appopriate names.

What is the controlDict dictionary?
===================================

Controls overall behaviour of run, timestep and saving behaviour

.. code-block:: c

    startFrom startTime;        // simulation starts at t=0
    startTime 0;
    
    stopAt endTime;             // simulation ends at t=10
    endTime 10;
    
    deltaT 0.005;               // timestep = 0.005s
    
    writeControl timeStep;      // write every 100 timesteps (uncompressed ASCII format)
    writeInterval 100;

**latestTime = continue after you have started simulation for a few seconds**    
    
What is the polyMesh dictionary?
================================

Stores the mesh structure containing (from blockMeshDict):

* Boundary - list of boundary types 
* Points - list of the mesh vertices
* Faces - list of which vertices make up which faces
* Owner - list of which cells own which faces
* Neighbour - list of which face has which neighbour cell

Located in ``case/constant/polyMesh``

At start of run OpenFOAM reads this information, checks it and constructs a mesh

What is the blockMeshDict dictionary?
=====================================

* Lists boundary patches
* Used by ``$ blockMesh`` utility to generate mesh (creates polyMesh dictionary)
* Located in ``case/system``

How are vertices created in blockMeshDict dictionary?
=====================================================

Right hand rule starting bottom left

.. code-block:: c 

    vertices
    (
        (0 0 0)             // node 0
        (1 0 0)             // node 1
        (1 1 0)             // node 2
        (0 1 0)             // node 3
        (0 0 0.1)           // node 4
        (1 0 0.1)           // node 5
        (1 1 0.1)           // node 6
        (0 1 0.1)           // node 7
    );

How are blocks created in blockMeshDict dictionary?
===================================================

.. code-block:: c

    blocks
    (
        hex (0 1 2 3 4 5 6 7) (20 20 1) simpleGrading (1 1 1)   // vertex numbers, cells (x y z), expansion ratios (x y z)
    );

How are patches created in blockMeshDict dictionary?
====================================================

* Patches = boundary type, name and location
* Points must be connected in sequence  - right hand rule 


.. code-block:: c

    boundary
    (
        movingWall
        {
            type wall;
            faces
            (
                (3 7 6 2)
            );
        }
        fixedWalls
        {
            type wall;
            faces
            (
                (0 4 7 3)
                (2 6 5 1)
                (1 5 4 0)
            );
        }
        frontAndBack
        {
            type empty;     // no calculations in z-direction, OpenFOAM is a 3D code, 2D is done with single thickness mesh in z-dir
            faces
            (
                (0 3 2 1)
                (4 5 6 7)
            );
        }
    );

How is pressure initialised in the ``0/p`` dictionary?
======================================================

.. code-block:: c

    dimensions      [0 2 -2 0 0 0 0];       // all fields have dimensions

    internalField   uniform 0;              // uniform internal field

    boundaryField                           // boundary conditions: patch name { type type info }
    {                                       // ordering same as blockMeshDict
        movingWall
        {
            type            zeroGradient;
        }

        fixedWalls
        {
            type            zeroGradient;
        }

        frontAndBack
        {
            type            empty;
        }
    }

**Pressure in OpenFOAM is pressure divided by density, units $m^2 s^{-2}$**

How is velocity initialised in the ``0/U`` dictionary?
======================================================

.. code-block:: c

    dimensions      [0 1 -1 0 0 0 0];           // m/s units

    internalField   uniform (0 0 0);            // u, v, w = 0

    boundaryField
    {
        movingWall
        {
            type            fixedValue;         // constant velocity = 1
            value           uniform (1 0 0);
        }

        fixedWalls
        {
            type            fixedValue;         // no slip
            value           uniform (0 0 0);
        }

        frontAndBack                            // no calculations
        {
            type            empty;
        }
    }

    
How are discretisation schemes specified in the fvScheme dictionary?
====================================================================

For example ``div(phi,q)`` - this requires an entry in fvSchemes

.. code-block:: c

    divSchemes
    (
        default none;
        div(phi,q) Gauss upwind;    // discretisation method is always Gauss
                                    // interpolation method is upwind (can be another method)
    )

Which interpolation methods are avaliable in the fvScheme dictionary?    
=====================================================================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Scheme
     - Meaning
   * - ::

           upwind

     - 1st order upwind
   * - ::

           linearUpwind

     - 2nd order correction to upwind  
   * - ::

           linearUpwindV

     - Improved handling for vectors
   * - ::

           linear

     - Central differencing (2nd order)
   * - ::

           SFCD

     - Self filtered central differencing     
   * - ::

           vanLeer

     - van Leer limited CD     

    
How are temporal schemes specified in the fvScheme dictionary?
==============================================================

.. code-block:: c

    ddtSchemes
    (
        default Euler;
    )


Which temporal schemes are avaliable in the fvScheme dictionary?    
================================================================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Scheme
     - Meaning
   * - ::

           steadyState

     - Steady state
   * - ::

           Euler

     - Euler
   * - ::

           CrankNicholson

     - Crank-Nicholson
   * - ::

           backward

     - Backward differencing
     

Which solvers are avaliable in the fvSolution dictionary, in the solvers subdictionary?    
=======================================================================================

Specify which solver to use for which equation (U, p) etc.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Meaning
   * - ::

           smoothSolver

     - Solver with smoothing (for speed)
   * - ::

           PCG
           

     - Preconditioned conjugate gradient (for stability) 
   * - ::

           PBiCG

     - Preconditioned biconjugate gradient (for stability)
   * - ::

           GAMG

     - Algebraic multigrid (for speed)
    
The solvers subdictionary also specifies:

* The tolerance and relative tolerance (relative tolerance = 0 for transient cases)
* Smoother and preconditioners

Which how is PISO specified in the fvSolution dictionary, in the PISO subdictionary?    
====================================================================================

Defines:

* ``nCorrectors`` - number of correction iterations (usually 2)
* ``nNonOrthogonalCorrectors`` - if mesh is badly non-orthogonal
* ``pRefCell`` and ``pRefValue`` - for reference values, if no boundaries have pressure

Which how is SIMPLE specified in the fvSolution dictionary, in the SIMPLE subdictionary?    
========================================================================================

Defines:

* ``nNonOrthogonalCorrectors`` - if mesh is badly non-orthogonal (not normally neccessary)

How are relaxation factors specified in the fvSolution dictionary, in the relxationFactors subdictionary?    
=========================================================================================================

Defines:

* Underrelaxation factors for equations (typically 0.3 to 0.9)

Lower relaxation factors means:

* Slower convergence
* More stable

How are iterations introduced in the controlDict dictionary?
============================================================

* set ``deltaT = 1`` in the controlDict



Utilities
---------

What are utilities?
===================

Utilities are used for:

* Mesh generation, conversion, manipulation
* Pre/post processing
* Data conversion


What are the mesh generation utilities?
=======================================

::

    $ blockMesh         // block structured mesh
    $ snappyHexMesh     // hexa mesh truncated at boundaries defined by STL files
    $ fluentMeshToFoam  // .msh
    $ gambitToFoam      // .neu
    $ startToFoam       // .neu
    $ ideasToFoam       // .ans
    $ cfxToFoam         // .geo
    
What are the post-processing utilities?
=======================================

::

    $ sample    // sample data along a line between two points
    $ flowType  // calculates a flow parameter
    $ yPlusRAS  // calculates y plus on walls
    $ yPlusLES  // calculates y plus on walls
    $ foamCalc  // takes a number of arguments for general field operations e.g. foamCalc mag U

What are functionObjects?
=========================

Evalautes quantities during run, by adding functions to controlDict e.g.

* Evaulate forces
* Manipulate fields
* Sample data
* Visualise flows
* Control execution
    
e.g.

.. code-block:: c

    forces
    {
        type forceCoeffs;
    }

    
Solvers
-------

What are solvers?
=================

Solve computational continuum mechanics problems e.g. turbulent flow, stress analysis


    
How are solvers invoked?
========================

e.g. to run the icoFoam solver:

::

    $ icoFoam

What are the simple solvers?
============================
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           laplacianFoam

     - Laplace Equation
   * - ::

           potentialFoam

     - Potential flow
   * - ::

           scalarTransportFoam
           
     - Tranpsort equation for a given velocity field
   * - ::

           icoFoam
           
     - Transient, incompressible, laminar, Newtonian flow (based on PISO)
   * - ::

           nonNewtonianIcoFoam
           
     - Transient, incompressible, laminar, non-Newtonian flow (based on PISO)       
   * - ::

           sonicFoam
           
     - Transonic/supersonic transient gas flow      
     
What are the turbulent solvers?
===============================
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           simpleFoam

     - Steady state (SIMPLE) solver for turbulent flows
   * - ::

           pisoFoam

     - Transient (PISO) solver for turbulent flows
   * - ::

           pimpleFoam
           
     - Large time step transient solver using merged PISO-SIMPLE algorithm (very stable with large timestep)
   * - ::

           pimpleDymFoam
           
     - Same as pimpleFoam but with mesh motion
   * - ::

           boundaryFoam
           
     - Steady state 1D turbulent flow used to generate boundary layer conditions at inlet      
   * - ::

           channelFoam
           
     - LES code for cyclic channels 
   * - ::

           porousSimpleFoam
           
     - Turbulent flow in a porous medium 

(The laminar option is avaliable in turbulent codes)

What are the multi-phase solvers?
=================================
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           bubbleFoam
           twoPhaseEulerFoam

     - Solvers for dispersed phase flow, e.g. gas bubbles in liquid
   * - ::

           interFoam
           multiphaseInterFoam

     - Solves for 2 (or multiple) immiscible phases with interface capturing using VOF method, laminar flow case
   * - ::

           twoLiquidMixingFoam

     - Solver for two immiscible fluids
   * - ::

           settlingFoam

     - Settling of solid particles in liquid
          
What are the combustion solvers?
================================
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           XiFoam

     - Premixed/partially premixed combustion with RANS/LES turbulence modelling and Weller model
   * - ::

           engineFoam
           coldEngineFlow

     - Solver for IC engine calculations with/without combustion            
   * - ::

           dieselFoam
           dieselEngineFoam

     - Diesel engine spray/combustion codes         
   * - ::

           reactingFoam
           rhoReactingFoam

     - Chemical reaction code    
     
What are the other solvers?
===========================
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           buoyant...Foam

     - Transient/s.s. solver comp. fluid, h.t. w/wo Boussinesq
   * - ::

           chtMultiRegionFoam

     - CHT between solid and fluid regions          
   * - ::

           dsmcFoam

     - Direct simulation Monte Carlo solver for multi-species flow         
   * - ::

           mdFoam

     - Molecular dynamics solver
   * - ::

           solidDisplacementFoam

     - Transient/s.s. solver for linear-elastic deformation
   * - ::

           financialFoam

     - Black-Scholes equations
    
    
    
What are the RANS solvers?
==========================
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           pisoFoam

     - Transient incompressible turbulent flow using PISO algorithm
   * - ::

           simpleFoam

     - Steady state incompressible turbulent flow using SIMPLE algorithm         
   * - ::

           spalartAllmaras
           kEpsilon
           realizableKE
           RNGkEpsilon
           NonlinearKEShih
           LienCubicKE
           kOmegaSST

     - 1 and 2 equation models
   * - ::

           LienCubicKEKLowRE
           LienLeschZinerLowRE
           LaunderSharmaKE
           LamBremhorstKE
           QZeta

     - Low Re models    
   * - ::

           LRR
           LaunderGibsonRSTM

     - Reynolds Stress models               

How are turbulence models invoked?
==================================

By a pointer, can switch between models with runtime selection.

::

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
        + fvm::div(phi,U)
        + turbulence -> divDevReff(U)
    );

    
How can y-plus be checked?
==========================

The following post-processing utility

::

    $ checkYPlus    
    
    
     
What are the mesh motion solvers?
=================================

Solvers include Dym in the name
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           dynamicFvMesh

     - Generic mesh motion
   * - ::

           pimpleDymFoam

     - Incompressible turbulent flow solver for moving mesh        
   * - ::

           rhoCentralDymFoam

     - Transonic (density-based) flow solver with moving mesh      
   * - ::

           sonicDymFoam

     - Sonic (pressure-based)
   * - ::

           interDymFoam

     - VOF multiphase solver
     

How is mesh motion specified?
=============================

In constant/dynamicMeshDict e.g.

* displacmentLaplacian
* velocityLaplacian

Different solver options use different fields, e.g.

* pointMotionU
* velocityLaplacian

     
What are the fluid-structure interaction solvers?
=================================================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           solidDisplacementFoam
           solidEquilibriumDisplacementFoam

     - Couping between fluid and solid

     
What are the free surface solvers?
==================================
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           interFoam

     - Basic VOF code (2 fluids) turbulence model can be RANS or LES
   * - ::

           compressibleInterFoam

     - Compressible flow version of interFoam         
   * - ::

           interDymFoam

     - Moving mesh version of interFoam       
   * - ::

           multiPhaseInterFoam

     - Solve for arbitary number of immiscible fluids
   * - ::

           cavitatingFoam

     - Allows for cavitation


What are the compressible flow solvers?
=======================================
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Solver
     - Purpose
   * - ::

           rhoSimpleFoam
           rhoPimpleFoam
           rhoPisoFoam

     - Compressible versions of standard OpenFOAM codes (pressure-based)
   * - ::

           rhoCentralFoam

     - MUSCL scheme of Kuranov and Tadmov         
   * - ::

           sonicFoam

     - Density based scheme for high speed compressible flows (transient solution)        
   * - ::

           sonicLiquidFoam

     - Same but for compressible liquids

All codes capable of laminar and turbulent simulation                  
  
    
Example - Boussineq Approximation
---------------------------------

For buoyancy-driven flow we often make use of the Boussinesq approximation : air modelled as incompressible with a body force proportional to ∆θ. Can we implement this into icoFoam?
Need to solve standard heat conduction equation:

.. math::
    {{\partial \theta} \over {\partial t}} + \nabla \cdot (u \theta) = {\kappa \over {\rho_0 C_v}} \nabla^2 \theta
    
and alter the momentum equation:

.. math::
    {{\partial u} \over {\partial t}} + \nabla \cdot(u u) = -\nabla p + \nu \nabla^2 u - \beta g (\theta_0 - \theta)
 
to accommodate this.

1. Copy ``icoFoam`` solver to user directory
    * ``cd $FOAM_RUN/../applications``
    * ``cp -r $WM_PROJECT_DIR/applications/solvers/incompressible/icoFoam .``
2. Compile ``icoFoam`` in user directory
    * ``cd icoFoam/Make``
    * edit ``files`` to change ``$FOAM_APPBIN`` to ``$FOAM_USER_APPBIN``
    * ``cd ..``
    * ``wmake``
3. Understand standard icoFoam (OpenFOAM 2.4.0)

.. code-block:: c
   
        #include "fvCFD.H" // This file brings in the most fundamental tools for performing any finite volume calculation and includes a bunch of other files.

        int main(int argc, char *argv[])
        {
            #include "setRootCase.H" // Checks folder structure of the case
            #include "createTime.H" // Checks runtime according to the controlDict and initiates time variables
            #include "createMesh.H" // Defines mesh in the domain
            #include "createFields.H" // Creates fields (e.g. u, p, T) according to createFields.H 
            #include "initContinuityErrs.H" // Declare and initialise the cumulative continuity error

            Info<< "\nStarting time loop\n" << endl;

            while (runTime.loop()) // Start time loop
            {
                Info<< "Time = " << runTime.timeName() << nl << endl; // Print the current time

                #include "readPISOControls.H" // Read control parameters from fvSchemes
                #include "CourantNo.H" // Calculates and outputs the mean and maximum Courant Numbers

                // Set up the linear algebra for the momentum equation.  
                // The flux of U, phi, is treated explicity
                // using the last known value of U.

                fvVectorMatrix UEqn
                (
                    fvm::ddt(U)
                  + fvm::div(phi, U)
                  - fvm::laplacian(nu, U)
                );
                
                // Solve using the last known value of p on the RHS.  
                // This gives us a velocity field that is
                // not divergence free, but approximately satisfies momentum.  
                // See Eqn. 7.31 of Ferziger & Peric (predicted velocity)
                
                solve(UEqn == -fvc::grad(p));

                // --- PISO loop

                for (int corr=0; corr<nCorr; corr++)
                {
                
                    // From the last solution of velocity, extract the diag. term from the matrix and store the reciprocal
                    // Note that the matrix coefficients are functions of U due to the non-linearity of convection.

                    volScalarField rAU(1.0/UEqn.A());

                    // take a Jacobi pass and update U.  See Hrv Jasak's thesis eqn. 3.137 and Henrik Rusche's thesis, eqn. 2.43
                    // UEqn.H is the right-hand side of the UEqn minus the product of (the off-diagonal terms and U).
                    // Note that since the pressure gradient is not included in the UEqn. above, this gives us U without
                    // the pressure gradient.  Also note that UEqn.H() is a function of U.

                    volVectorField HbyA("HbyA", U);
                    
                    // Calculate the fluxes by dotting the interpolated velocity (to cell faces) with face normals
                    // The ddtPhiCorr term accounts for the divergence of the face velocity field by taking out the 
                    // difference between the interpolated velocity and the flux.
                    
                    HbyA = rAU*UEqn.H();
                    surfaceScalarField phiHbyA
                    (
                        "phiHbyA",
                        (fvc::interpolate(HbyA) & mesh.Sf())
                      + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                    );
                    
                    // Adjusts the inlet and outlet fluxes to obey continuity, which is necessary for creating a well-posed
                    // problem where a solution for pressure exists.
                    
                    adjustPhi(phiHbyA, U, p);

                    // Iteratively correct for non-orthogonality.  The non-orthogonal part of the Laplacian is calculated from the most recent
                    // solution for pressure, using a deferred-correction approach.

                    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                    {
                    
                        // Set up the pressure equation
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                        );

                        // In incompressible flow, only relative pressure matters.  Unless there is a pressure BC present,
                        // one cell's pressure can be set arbitrarily to produce a unique pressure solution

                        pEqn.setReference(pRefCell, pRefValue);
                        pEqn.solve();

                        // On the last non-orthogonality correction, correct the flux using the most up-to-date pressure
                        // The .flux method includes contributions from all implicit terms of the pEqn (the Laplacian)

                        if (nonOrth == nNonOrthCorr)
                        {
                            phi = phiHbyA - pEqn.flux();
                        }
                    } // end of non-orthogonality looping

                    #include "continuityErrs.H"

                    // Add pressure gradient to interior velocity and BC's.  
                    // Note that this pressure is not just a small correction to a previous pressure, 
                    // but is the entire pressure field.  
                    // Contrast this to the use of p' in Ferziger & Peric, Eqn. 7.37.

                    U = HbyA - rAU*fvc::grad(p);
                    U.correctBoundaryConditions();
                }

                runTime.write();

                Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                    << nl << endl;
            } // end of the time step loop

            Info<< "End\n" << endl;

            return 0;
            }
   
4. Modify icoFoam in the following way:
Open ``createFields.H`` and read in the various properties:

.. code-block:: c

    dimensionedScalar kappa
    (
        transportProperties.lookup ("kappa")
    );
    
(similar lines for rho0, Cv, theta0 and beta). Also worth introducing hCoeff:

.. code-block:: c

    dimensionedScalar hCoeff = kappa/( rho0 * Cv );
    
5. Introduce gravitational acceleration g; read in from the same dictionary, but is a dimensionedVector rather than a dimensionedScalar.

6. Create a temperature field theta as a volScalarField and read it in. This is very similar to the pressure field, so make a copy of this and modify accordingly.

7. Modify the momentum equation (UEqn) to add the term

.. code-block:: c

    + beta * g *( theta0 - theta )

8. Outside of the PISO loop – create and solve the temperature equation:

.. code-block:: c

    fvScalarMatrix tempEqn
    (
        fvm :: ddt (theta)
        + fvm :: div (phi, theta)
        - fvm :: laplacian (hCoeff, theta)
    );
    tempEqn.solve ();

9. Compile this using wmake (rename executable as boussinesqFoam) - folder name, icoFoam.C and files (icoFoam.C and icoFoam)


Example - Casson model
----------------------

The Casson model is a non-Newtonian viscosity model – used for chocolate, blood, etc. Can we implement in OpenFOAM? Stress-strain relation for a fluid

.. math::

    \tau = \rho \nu \dot{\gamma}

where
    
.. math::

    \dot{\gamma} = {1 \over 2}(\nabla \vec{u} + \nabla \vec{u}^T)

is rate of strain tensor

:math:`\nu` = const is a Newtonian fluid. :math:`\nu = \nu(\cdot{\gamma})` is non-Newtonian.

Casson model:

.. math::

    \nu(J_2) = {{[(\eta^2 J_2)^{1/4} + \sqrt{\tau_y / 2}]^2} \over {\rho \sqrt{J_2}}}

where

.. math::

    J_2 = ||{1 \over 2}(\nabla \vec{u} + \nabla \vec{u}^T)||^2

How to make -builtin create .foam extension files?
--------------------------------------------------

In:

::

    /OpenFOAM/OpenFOAM-2.4.0/bin/paraFoam 

Change these lines:

Line 71: 

::

    extension=foam
    
Line 181:

::
 
    extension=foam
    
Line 241:

::
 
    builtin | foam)
    
In:

::

    /OpenFOAM/OpenFOAM-v1612+/bin/paraFoam 

Change these lines:

Line 72: 

::

    extension=foam
    
Line 187:

::
 
    extension=foam


How to load new results into Paraview?
--------------------------------------

The first command creates .foam
The second comment opens Paraview

::

    $ paraFoam -touch
    $ paraFoam &
    
How to save a state?
--------------------

File > Save State

Filename.pvsm     
    
How to load a state?
--------------------

::

    $ paraview &
    
File > Load State  
    
Filename.pvsm    
    
        
Example - Mesh Refinement (OpenFOAM 2.4.0)
------------------------------------------

Jozsef Nagy 1, 2

::

    $ two
    $ run
    $ cp -r $FOAM_TUTORIALS/incompressible/icoFoam/elbow/ .
    $ mv elbow 001_elbow_tri
    $ cp -r 001_elbow_tri 002_elbow_quad
    $ cp -r 002_elbow_quad 003_elbow_quad_refined
    $ rm -rf 002_elbow_quad/elbow.msh
    $ rm -rf 003_elbow_quad_refined/elbow.msh
    
* Download elbow_quad.msh and copy to 002 and 003 directories

* Change endTime to 75s in controlDict for all cases

::

    $ cd 001_elbow_tri
    $ fluentMeshToFoam elbow.msh
    $ paraFoam &
    $ icoFoam
    $ paraFoam -touch
    $ paraFoam &
    
::

    $ cd ../002_elbow_quad    
    $ fluentMeshToFoam elbow_quad.msh
    $ icoFoam
    $ paraFoam -touch
    $ paraFoam &

::

    $ cd ../003_elbow_quad_refined 
    $ fluentMeshToFoam elbow_quad.msh
    $ refineMesh -overwrite
    
* Halve timestep deltaT in controlDict because spatial step has been halved
* Double writeInterval in controlDict    

::

    $ icoFoam
    $ paraFoam -touch
    $ paraFoam &

Create contour plots
    
* Load results from ``001_elbow_tri.OpenFOAM`` and ``002_elbow_quad.OpenFOAM``    
* Translate meshes - can use search to find translate
* Toggle legend on (first icon from left) 
* Edit color map (second icon from left)     
* Rescale to custom range (second icon down RHS) - change max to 4 - Rescale
* Toggle legend off and on
* Edit color legend properties (top left icon)
* Title = Velocity (m/s)
* Delete component title
* Change font to 10
* Apply - ok
* Can move legend so it's horizontal
* File > Save Screenshot > Ok > velocity.png

Create line plot

* Select the plot
* Filters > Data Analysis > Plot Over Line > Apply
* Move points to correspond with line of interest
* Use cursor to select the resulting chart (can change the variables that are plotted)
* File > Save Data > Enter filename
* If write all timesteps is not selected, then it will write only the last timestep.
* Open in Libre Office and plot

Example - Block Mesh (OpenFOAM 2.4.0)
-------------------------------------

Jozsef Nagy 3

Tab will fill in the rest of the folder.

::

    $ two
    $ run
    $ cp -r $FOAM_TUTORIALS/compressible/sonicFoam/laminar/forwardStep/ .
    $ mv forwardStep 004_forward_step

Connect to sftp server.

Run through of blockMeshDict

.. code-block:: c

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.4.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      blockMeshDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    convertToMeters 1;          // 1= values given in meters, 0.001 = millimeters

    vertices                    // it doesn't matter about order here
    (
        (0 0 -0.05)
        (0.6 0 -0.05)
        (0 0.2 -0.05)
        (0.6 0.2 -0.05)
        (3 0.2 -0.05)
        (0 1 -0.05)
        (0.6 1 -0.05)
        (3 1 -0.05)
        (0 0 0.05)
        (0.6 0 0.05)
        (0 0.2 0.05)
        (0.6 0.2 0.05)
        (3 0.2 0.05)
        (0 1 0.05)
        (0.6 1 0.05)
        (3 1 0.05)
    );

    blocks                      // order must be right hand rule (x y z) grading
    (
        hex (0 1 3 2 8 9 11 10) (25 10 1) simpleGrading (1 1 1)
        hex (2 3 6 5 10 11 14 13) (25 40 1) simpleGrading (1 1 1)
        hex (3 4 7 6 11 12 15 14) (100 40 1) simpleGrading (1 1 1)
    );

    edges
    (
    );

    boundary                    // right hand rule again - normal vector must point outwards
    (
        inlet
        {
            type patch;
            faces
            (
                (0 8 10 2)
                (2 10 13 5)
            );
        }
        outlet
        {
            type patch;
            faces
            (
                (4 7 15 12)
            );
        }
        bottom
        {
            type symmetryPlane;
            faces
            (
                (0 1 9 8)
            );
        }
        top
        {
            type symmetryPlane;
            faces
            (
                (5 13 14 6)
                (6 14 15 7)
            );
        }
        obstacle
        {
            type patch;
            faces
            (
                (1 3 11 9)
                (3 4 12 11)
            );
        }
    );

    mergePatchPairs
    (
    );

    // ************************************************************************* //

Run blockMesh

::

    $ blockMesh


Check boundary file in constant/polyMesh

.. code-block:: c

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.4.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       polyBoundaryMesh;
        location    "constant/polyMesh";
        object      boundary;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    6
    (
        inlet
        {
            type            patch;
            nFaces          50;
            startFace       10325;
        }
        outlet
        {
            type            patch;
            nFaces          40;
            startFace       10375;
        }
        bottom
        {
            type            symmetryPlane;
            inGroups        1(symmetryPlane);
            nFaces          25;
            startFace       10415;
        }
        top
        {
            type            symmetryPlane;
            inGroups        1(symmetryPlane);
            nFaces          125;
            startFace       10440;
        }
        obstacle
        {
            type            patch;
            nFaces          110;
            startFace       10565;
        }
        defaultFaces                        // unspecified front and back planes
        {
            type            empty;
            inGroups        1(empty);
            nFaces          10500;
            startFace       10675;
        }
    )

    // ************************************************************************* //

Open mesh in paraView

::

    $ paraFoam -touch
    $ paraFoam &

Change the mesh to wireframe, black lines and white background:

* Properties > Display > Representation > Wireframe
* Properties > Coloring > Solid Color > Edit > Black
* Edit > View Settings > General > Solid Color = White > Apply

Change default color space to RGB:

* Edit color map > Choose preset > Blue to Red Rainbow > Save as default

Pressure ICs and BCs

.. code-block:: c

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.4.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      p;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [1 -1 -2 0 0 0 0];

    internalField   uniform 1;                  // 1 Pa

    boundaryField
    {
        inlet
        {
            type            fixedValue;
            value           uniform 1;          // 1 Pa
        }

        outlet
        {
            type            waveTransmissive;
            field           p;
            phi             phi;
            rho             rho;
            psi             thermo:psi;
            gamma           1.4;
            fieldInf        1;
            lInf            3;
            value           uniform 1;
        }

        bottom
        {
            type            symmetryPlane;
        }

        top
        {
            type            symmetryPlane;
        }

        obstacle
        {
            type            zeroGradient;       // Neumann boundary condition (derivative)
        }

        defaultFaces
        {
            type            empty;
        }
    }

    // ************************************************************************* //

Temperature ICs and BCs

.. code-block:: c   
    
    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.4.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      T;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 0 0 1 0 0 0];

    internalField   uniform 1;              // 1K(!)

    boundaryField
    {
        inlet
        {
            type            fixedValue;
            value           uniform 1;      // 1K(!)
        }

        outlet
        {
            type            zeroGradient;   // Neumann
        }

        bottom
        {
            type            symmetryPlane;
        }

        top
        {
            type            symmetryPlane;
        }

        obstacle
        {
            type            zeroGradient;
        }

        defaultFaces
        {
            type            empty;
        }
    }

    // ************************************************************************* //

Velocity ICs and BCs

.. code-block:: c 

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.4.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volVectorField;
        object      U;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 1 -1 0 0 0 0];

    internalField   uniform (3 0 0);            // 3 m/s

    boundaryField
    {
        inlet
        {
            type            fixedValue;
            value           uniform (3 0 0);    // 3 m/s
        }

        outlet
        {
            type            zeroGradient;
        }

        bottom
        {
            type            symmetryPlane;
        }

        top
        {
            type            symmetryPlane;
        }

        obstacle
        {
            type            fixedValue;         // No slip
            value           uniform (0 0 0);
        }

        defaultFaces
        {
            type            empty;
        }
    }

    // ************************************************************************* //

Turbulence Properties - laminar

.. code-block:: c

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.4.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      turbulenceProperties;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    simulationType  laminar;


    // ************************************************************************* //

Themophyscial properties

.. code-block:: c

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.4.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      thermophysicalProperties;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    thermoType
    {
        type            hePsiThermo;
        mixture         pureMixture;
        transport       const;
        thermo          hConst;
        equationOfState perfectGas;
        specie          specie;
        energy          sensibleInternalEnergy;
    }

    // Note: these are the properties for a "normalised" inviscid gas
    //       for which the speed of sound is 1 m/s at a temperature of 1K
    //       and gamma = 7/5
    mixture
    {
        specie
        {
            nMoles          1;
            molWeight       11640.3;
        }
        thermodynamics
        {
            Cp              2.5;        // Specific heat
            Hf              0;
        }
        transport
        {
            mu              0;
            Pr              1;          // Prandtl number
        }
    }


    // ************************************************************************* //


Write Interval at 0.5 seconds specified using ``runTime``

.. code-block:: c

    /*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  2.4.0                                 |
    |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      controlDict;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    application     sonicFoam;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         10;

    deltaT          0.002;

    writeControl    runTime;

    writeInterval   0.5;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;


    // ************************************************************************* //

Run simulation

::

    $ sonicFoam


* Load paraView:

::

    $ paraFoam -builtin &    


Change the colors to black and white:

Edit > Settings > Colors > Palette = Print > Apply > Ok     
    
Can refresh if ParaView was already open:

* Properties > Refresh    
* Toggle: Surface, Wireframe
    
Change grading:

.. code-block:: c

    blocks                      // order must be right hand rule (x y z) grading
    (
        hex (0 1 3 2 8 9 11 10) (25 10 1) simpleGrading (0.5 1 1)
        hex (2 3 6 5 10 11 14 13) (25 40 1) simpleGrading (0.5 1 1)
        hex (3 4 7 6 11 12 15 14) (100 40 1) simpleGrading (1 1 1)
    );
    
Remove all folders except ICS:

::

    $ rm -rf 0.* [1-9]*
    
Overwrite files in constant/polyMesh

:: 

    $ blockMesh        
    
Reload mesh:
    
* Properties > Refresh    
* Toggle: Surface, Wireframe  

The expansion ratio is:

.. math::

    \text{Ratio} = {\text{End spacing} \over \text{Start spacing}}
    
With the directions given by the Cartesian coordinates
    
Can change to edge spacing:

.. code-block:: c

    blocks                      // order must be right hand rule (x y z) grading
    (
        hex (0 1 3 2 8 9 11 10) (50 20 1) edgeGrading (1 1 1 1 2 2 2 2 1 1 1 1)
        hex (2 3 6 5 10 11 14 13) (25 40 1) simpleGrading (1 1 1)
        hex (3 4 7 6 11 12 15 14) (100 40 1) simpleGrading (1 1 1)
    );  
    
Overwrite files in constant/polyMesh

:: 

    $ blockMesh        
    
Reload mesh:
    
* Properties > Refresh    
* Toggle: Surface, Wireframe     
    
Example - Grid Convergence (OpenFOAM 2.4.0)
-------------------------------------------

Jozsef Nagy 4

Tab will fill in the rest of the folder.

::

    $ two
    $ run
    $ cp -r $FOAM_TUTORIALS/compressible/sonicFoam/laminar/shockTube/ .
    $ mv shockTube 005_shock_tube
    $ cd 005_shock_tube
    $ gedit ./constant/polyMesh/blockMeshDict
    $ cd ..
    
Change the number of x cells from 1000 to 100 and save.

Change the name and copy the folders for 1000 and 10000 cells

::

    $ mv 005_shock_tube 005_shock_tube_100
    $ cp -r 005_shock_tube_100 006_shock_tube_1000
    $ cp -r 005_shock_tube_100 007_shock_tube_10000
    
Create the mesh for the 100 cell case

::

    $ cd 005_shock_tube_100
    $ blockMesh
  
Open case in Paraview

::

    $ paraFoam -builtin &
    
Remove region in setFieldsDict in order to set all the fields to be the same

::

    $ nano ./system/setFieldsDict
    
.. code-block:: c    
    
    defaultFieldValues ( volVectorFieldValue U ( 0 0 0 ) volScalarFieldValue T 348.432 volScalarFieldValue p 100000 );

    regions         ( 
    // boxToCell { box ( 0 -1 -1 ) ( 5 1 1 ) ; fieldValues ( volScalarFieldValue T 278.746 volScalarFieldValue p 10000 ) ; } 
    );       
    
Now set the fields - this sets the fields of 0/p, 0/U, 0/magU and 0/T to uniform values

::

    $ setFields
    
Now set a non-uniform field from centre to end in x and all y and all z (RHS of geometry)

::

    $ nano ./system/setFieldsDict


.. code-block:: c    
    
    defaultFieldValues ( volVectorFieldValue U ( 0 0 0 ) volScalarFieldValue T 348.432 volScalarFieldValue p 100000 );

    regions         ( 
    boxToCell { box ( 0 -1 -1 ) ( 5 1 1 ) ; fieldValues ( volScalarFieldValue T 278.746 volScalarFieldValue p 10000 ) ; } 
    );     

ctrl + o
ctrl + x

This sets the fields of 0/p, 0/U, 0/magU and 0/T to non-uniform values 

::

    $ setFields

Check that the fields are correct:

::

    $ nano ./0/p

Run the simulation

::

    $ sonicFoam

Refresh the results in ParaView

* Properties > Refresh    
* Toggle: Surface, Wireframe

Plot the values as a line chart:

* Filters > Search > plot over line > Enter
* Apply
* Change y and z values to 0 for centreline > Apply
* Top right order of chart will change from fullscreen to split screen and back
* Shows that grid density is too coarse
* File > Save State > name.pvsm

(You have to run blockMesh on the grid before loading a state paraview, as foam.foam is empty - but you can load foam.foam from stratch)

Now simulate different grids

::

    $ cd ../006_shock_tube_1000

Change resolution to 1000 in blockMeshDict

::

    $ nano ./constant/polyMesh/blockMeshDict

ctrl + o
ctrl + x

Create mesh and set fields

::

    $ blockMesh
    $ setFields

Check that the fields are correct:

::

    $ nano ./0/p

Run:

::

    $ sonicFoam


Load in Paraview

::

    $ paraFoam -builtin &

Compare with coarse case:

File > Open > X.foam from coarse directory

Change resolution to 10000 in blockMeshDict

::

    $ cd ../007_shock_tube_10000

::

    $ nano ./constant/polyMesh/blockMeshDict

ctrl + o
ctrl + x

Create mesh and set fields

::

    $ blockMesh
    $ setFields

Check that the fields are correct:

::

    $ nano ./0/p

Run:

::

    $ sonicFoam

Simulation will diverge - due to Courant number being too high (C > 1). Reducing cell spacing by a factor of 10 increases Courant number by a factor of 10 from previous value (0.28)

::

    Courant Number mean: 0.00264441 max: 2.84091

Reduce timestep by a factor of 10

::

    $ nano ./system/controlDict
	
Run:

::

    $ sonicFoam

Compare with coarse case:

File > Open > X.foam from coarse directory

Translate the grid by 2.1

Use sampleDict to sample 10,000 points on 10,000 cell grid.

::

	$ nano ./system/sampleDict

.. code-block:: c 

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	interpolationScheme cellPoint;

	setFormat       raw;

	sets
	(
	    data
	    {
		type    uniform;
		axis    x;
		start   (-4.995 0 0);
		end     (4.995 0 0);
		nPoints 10000;
	    }
	);

	fields          (T magU p);

	// ************************************************************************* //

To list all folders (should not see magU in all folders):

::

	$ tree

Calculate the magnitude of the velocity (magU should now appear):

::

	$ foamCalc mag U
	$ tree

Now sample the data:

::

	$ sample

Check if sample has created: 

::

	$ tree

Do the same for 1000 (10,000 points):

::

	$ cd ../006_shock_tube_1000
	$ nano ./system/sampleDict 

Magnitude of velocity:

::

	$ foamCalc mag U

Sample the velocity:

::

	$ sample

Do the same for 100 (10,000 points):

::

	$ cd ../005_shock_tube_100
	$ nano ./system/sampleDict 

Magnitude of velocity:

::

	$ foamCalc mag U

Sample the velocity:

::

	$ sample

Plot three mesh densities in one chart:

.. code-block:: python 

	import numpy as np
	import matplotlib.pyplot as plt
	import math as mt

	data10000 = np.genfromtxt('postProcessing/sets/0.007/data_T_magU_p.xy', delimiter=' ', skip_header=1,
		                 names=['x','T','magU','p'])
	x = data10000['x']
	T = data10000['T']
	magU = data10000['magU']
	p = data10000['p']

	data1000 = np.genfromtxt('../006_shock_tube_1000/postProcessing/sets/0.007/data_T_magU_p.xy', delimiter=' ', skip_header=1,
		                 names=['x','T','magU','p'])
	x = data1000['x']
	T = data1000['T']
	magU = data1000['magU']
	p = data1000['p']

	data100 = np.genfromtxt('../005_shock_tube_100/postProcessing/sets/0.007/data_T_magU_p.xy', delimiter=' ', skip_header=1,
		                 names=['x','T','magU','p'])
	x = data100['x']
	T = data100['T']
	magU = data100['magU']
	p = data100['p']

	f = plt.figure(1)

	ax1 = plt.gca()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	line1=ax1.plot(data10000['x'], data10000['p'], color='b', label='10000 cells')
	line2=ax1.plot(data1000['x'], data1000['p'], color='g', label='1000 cells')
	line3=ax1.plot(data100['x'], data100['p'], color='r', label='100 cells')

	# added these three lines
	lines = line1+line2+line3
	labels = [l.get_label() for l in lines]
	legend= ax1.legend(lines, labels, loc=0, fontsize='medium')

	ax1.set_xlabel(r'Distance, x \textit{(m)}')
	ax1.set_ylabel(r'Pressure \textit{(Pa)}')
	ax1.set_ylim([0,120000])
	plt.savefig('shock_tube_pressure')

	g = plt.figure(2)

	ax1 = plt.gca()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	line1=ax1.plot(data10000['x'], data10000['T'], color='b', label='10000 cells')
	line2=ax1.plot(data1000['x'], data1000['T'], color='g', label='1000 cells')
	line3=ax1.plot(data100['x'], data100['T'], color='r', label='100 cells')

	lines = line1+line2+line3
	labels = [l.get_label() for l in lines]
	legend= ax1.legend(lines, labels, loc=0, fontsize='medium')

	ax1.set_xlabel(r'Distance, x \textit{(m)}')
	ax1.set_ylabel(r'Temperature \textit{(K)}')
	ax1.set_ylim([200,500])
	plt.savefig('shock_tube_temperature')

	h = plt.figure(3)

	ax1 = plt.gca()
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	line1=ax1.plot(data10000['x'], data10000['magU'], color='b', label='10000 cells')
	line2=ax1.plot(data1000['x'], data1000['magU'], color='g', label='1000 cells')
	line3=ax1.plot(data100['x'], data100['magU'], color='r', label='100 cells')

	# added these three lines
	lines = line1+line2+line3
	labels = [l.get_label() for l in lines]
	legend= ax1.legend(lines, labels, loc=0, fontsize='medium')

	ax1.set_xlabel(r'Distance, x \textit{(m)}')
	ax1.set_ylabel(r'Velocity \textit{(m/s)}')
	ax1.set_ylim([-50,350])
	plt.savefig('shock_tube_velocity')

	plt.show()

Plot pressure, velocity and temperature:

.. code-block:: python

	import numpy as np
	import matplotlib.pyplot as plt
	import math as mt

	data = np.genfromtxt('postProcessing/sets/0.007/data_T_magU_p.xy', delimiter=' ', skip_header=1,
		                 names=['x','T','magU','p'])
	x = data['x']
	T = data['T']
	magU = data['magU']
	p = data['p']

	ax1 = plt.gca()
	ax2 = ax1.twinx()

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	line1=ax1.plot(data['x'], data['p'], color='b', label='Pressure (Pa)')
	line2=ax2.plot(data['x'], data['T'], color='r', label='Temperature (K)')
	line3=ax2.plot(data['x'], data['magU'], color='g', label='Velocity Magnitude (m/s)')

	lines = line1+line2+line3
	labels = [l.get_label() for l in lines]
	legend= ax1.legend(lines, labels, loc=0, fontsize='medium')

	ax1.set_xlabel(r'Distance, x \textit{(m)}')
	ax1.set_ylabel(r'Pressure \textit{(-)}')
	ax2.set_ylabel(r'Temperature or Velocity Magnitude \textit{(-)}')
	ax1.set_ylim([0,120000])
	ax2.set_ylim([-100,600])
	plt.savefig('005_shock_tube_100')
	plt.show()

Example - Transport Equation (OpenFOAM 2.4.0)
---------------------------------------------

Jozsef Nagy 5

Copy transport properties from pitz daily.

::

    $ two
    $ run
    $ cp -r $FOAM_TUTORIALS/compressible/sonicFoam/laminar/shockTube/ .
    $ cp -r $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/ .
	$ mv shockTube 008_tranport_base_case
	$ mv pitzDaily 009_pitz_daily
	$ cd 008_tranport_base_case
	$ rm ./0/magU ./0/p
	$ rm ./constant/thermophysicalProperties ./constant/turbulenceProperties
	$ cp -r ../009_pitz_daily/constant/transportProperties ./constant

Copy fvSchemes, fvSolution and controlDict from pitz daily

::

	$ rm ./system/controlDict ./system/fvSchemes ./system/fvSolution
	$ cp ../009_pitz_daily/system/controlDict ./system
	$ cp ../009_pitz_daily/system/fvSchemes ./system
	$ cp ../009_pitz_daily/system/fvSolution ./system

Change controlDict for an endTime = 5, writeControl = runTime, writeInterval=1

::

	$ nano ./system/controlDict 

Remove sampleDict

::

	$ rm ./system/sampleDict

Change setFieldsDict so U is 0, T is 0 by default and there is a region where T=1

::

	$ nano ./system/setFieldsDict

::

	defaultFieldValues ( volVectorFieldValue U ( 0 0 0 ) volScalarFieldValue T 0 );
	regions         ( boxToCell { box ( -0.5 -1 -1 ) ( 0.5 1 1 ) ; fieldValues ( volScalarFieldValue T 1 ) ; } );

Create the mesh and set the fields

::

	$ blockMesh
	$ setFields

Check U and T

::

	$ nano ./0/U
	$ nano ./0/T

Open Paraview and check temperature is initialised correctly

::

    $ paraFoam -touch
    $ paraFoam -builtin &

Make copies of base case

::

	$ cp -r 008_tranport_base_case 009_transport_base_case
	$ cp -r 008_tranport_base_case 010_transport_base_case
	$ cp -r 008_tranport_base_case 011_transport_base_case
	$ cp -r 008_tranport_base_case 012_transport_base_case
	$ cp -r 008_tranport_base_case 013_transport_base_case
	$ cp -r 008_tranport_base_case 014_transport_base_case
	$ cp -r 008_tranport_base_case 015_transport_base_case
	$ cp -r 008_tranport_base_case 016_transport_base_case
	$ cp -r 008_tranport_base_case 017_transport_base_case
	$ cp -r 008_tranport_base_case 018_transport_base_case
	$ cp -r 008_tranport_base_case 019_transport_base_case

Set the values like this

.. figure:: ../_images/nagy.png
   :scale: 75%
   :align: center

a) Run 008_transport_base_case

::

    $ scalarTransportFoam
    
    
b) Open in Paraview

::

    $ paraFoam -touch
    $ paraFoam -builtin
    
c) Save state

Do a), b) and c) for 009, 010, 011, 012, 013, 014, 015, 016, 017, 018 and 019_transport_base_case

Example - Discretisation (OpenFOAM 2.4.0)
-----------------------------------------

Orthogonal means the vector from the face is parallel to the vector joining the cell centres.

Non-orthogonal means the vector from the face is not parallel to the vector from the cell centres.

::

    $ two
    $ cp -r $FOAM_TUTORIALS/compressible/sonicFoam/laminar/shockTube/ .
    $ mv shockTube 019_discretisation_base_case
    $ cd 019_discretisation_base_case
    $ cd 0
    $ rm magU p
    $ nano T
    
Change sides to inletOutlet (from zeroGradient). Velocity is fixed to zero entering domain.

::

    sides
    {
        type        inletOutlet;
        inletValue  uniform 0;
        value       uniform 0;
    }

Copy pitzDaily transportProperties from scalarTransportFoam 

::

    $ cd ../constant
    $ rm thermophysicalProperties turbulenceProperties
    $ cp -r $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/constant/transportProperties .

Set diffusivity to zero - so that only convection will be used.

remove controlDict fvSchemes and fvSolution

::

    $ rm ./system/controlDict ./system/fvSchemes ./system/fvSolution
    
Copy files from pitzDaily for scalarTrasportFoam

::

    $ cp $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/system/controlDict ./system    
    $ cp $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/system/fvSchemes ./system 
    $ cp $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/system/fvSolution ./system 

Change endtime to 5 in controlDict, deltaT = 0.005 (for CFL<1 for Quick scheme), writeControl = runTime, writeInterval = 1;

::
    
    $ nano ./system/controlDict
    
Set the fields to a velocity of 1 in positive x and passive scalar to 1

::

    $ nano ./system/setFieldsDict
    
::

    defaultFieldValues ( volVectorFieldValue U ( 1 0 0 ) volScalarFieldValue T 0 );

    regions ( boxToCell { box ( -0.5 -1 -1 ) ( 0.5 1 1 ) ; fieldValues ( volScalarFieldValue T 1 ) ; } );

    
Create mesh

::

    $ blockMesh    
       
setFields

::

    $ setFields       
          
Create other cases:

::

    $ cd ..
    $ cp -r 019_discretisation_base_case 020_upwind
    $ cp -r 019_discretisation_base_case 021_linear
    $ cp -r 019_discretisation_base_case 022_linear_upwind
    $ cp -r 019_discretisation_base_case 023_quick
    $ cp -r 019_discretisation_base_case 024_cubic
	
Change the schemes

- Upwind for divergence term in transport equation

::

    $ cd 020_upwind
    $ nano ./system/fvSchemes

::

    divSchemes
    {
        ...
        div(phi,T) Gauss upwind;
    }

Start the simulation

::

    $ scalarTransportFoam

Create .foam file

::

    $ paraFoam -touch
    
Open Paraview

::

    $ paraFoam -builtin
    
Can rescale to range to observe diffusion

- Linear for divergence term in transport equation

::

    $ cd 021_linear
    $ nano ./system/fvSchemes

::

    divSchemes
    {
        ...
        div(phi,T) Gauss linear;
    }
    
Start simulation, create .foam file and open in paraview (from existing paraview)    

Linear scheme maintains 1, but now has -ve values - oscillations seen by rescaling from -0.01 to 0.01

The inletOutlet condition is one derived from mixed, which switches between zeroGradient when the fluid flows out of the domain at a patch face, and fixedValue, when the fluid is flowing into the domain. This was applied for temperature T.

The inletOutlet boundary condition is normally the same as zeroGradient, but it switches to fixedValue if the velocity vector next to the boundary aims inside the domain (backward flow). The value of that fixedValue is inletValue.

- Linear upwind for divergence term in transport equation

::

    $ cd 022_linear_upwind
    $ nano ./system/fvSchemes

The grad(T) is needed in the linear upwind differencing scheme - will use linear scheme defined for gradSchemes  
    
::

    gradSchemes
    {
        default Gauss linear;
    }

    divSchemes
    {
        ...
        div(phi,T) Gauss linearUpwind grad(T);
    }

Start simulation, create .foam file and open in paraview (from existing paraview)  

No more oscillations in T using linearUpwind, but still slightly negative value.

- QUICK for divergence term in transport equation

::

    $ cd 023_quick
    $ nano ./system/fvSchemes
  
::

    divSchemes
    {
        ...
        div(phi,T) Gauss QUICK;
    }

The minimum in QUICK is zero, QUICK and linearUpwind are comparable.    
    
- 4th order cubic for divergence term in transport equation

::

    $ cd 024_cubic
    $ nano ./system/fvSchemes
  
::

    divSchemes
    {
        ...
        div(phi,T) Gauss cubic;
    }   
    
Cubic has some oscillations

Change T boundary sides from inletOutlet to zeroGradient

Deleted old folders

::

    $ rm -rf *[1-9]
    
Re-run simulation

::
    
    $ scalarTransportFoam
    
Refresh in Paraview - shows large oscillations - wave is being reflected.

QUICK and linear upwind are the best.

Change Courant number to 1 for Linear Upwind deltaT = 0.01:

::

    nano system/controlDict
    
    
Deleted old folders

::

    $ rm -rf *[1-9]

Re-run

::

    $ scalarTransportFoam
    
Refresh in paraview - field is now more smeared, but it runs

Using Courant number 1 for QUICK - it doesn't run

Example - Blob (OpenFOAM 2.4.0)
-------------------------------

::

    $ two
    $ run
    $ cp -r $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/ .
    $ mv pitzDaily 025_blob
    $ cd 025_blob
    
Seeing what Pitz Daily actually does

::

    $ blockMesh
    $ scalarTransportFoam
    $ paraFoam -touch
    $ paraFoam -builtin
    
Realising that Pitz Daily has the wrong mesh and wrong boundary conditions

Copying the mesh from the lid driven cavity:

::

    $ cp -r /home/apr207/OpenFOAM/OpenFOAM-2.4.0/tutorials/incompressible/icoFoam/cavity/constant/polyMesh/blockMeshDict ./constant/polyMesh/

    
Changing the blockMeshDict to match the geometry:

::

    convertToMeters 1.0;

    vertices
    (
        (0 0 0)
        (10 0 0)
        (10 10 0)
        (0 10 0)
        (0 0 1)
        (10 0 1)
        (10 10 1)
        (0 10 1)
    );

    blocks
    (
        hex (0 1 2 3 4 5 6 7) (50 50 1) simpleGrading (1 1 1)
    );

    edges
    (
    );

    boundary
    (
        sides
        {
            type patch;
            faces
            (
                (3 7 6 2)
                (0 4 7 3)
                (2 6 5 1)
                (1 5 4 0)
            );
        }
        frontAndBack
        {
            type empty;
            faces
            (
                (0 3 2 1)
                (4 5 6 7)
            );
        }
    );

    mergePatchPairs
    (
    );

Possible setFields (tutorials/compressible/sonicFoam/laminar)

- boxToCell
- cylinderToCell
- sphereToCell

::

    cylinderToCell
    (
        p1 (5 5 0);
        p2 (5 5 1);
        radius 1;
        fieldValues
        (
            volScalarFieldValue T 600
        );
    )

    sphereToCell
    (
        centre (0 0 0);
        radius 1;
        fieldValues
        (
            volScalarFieldValue T 600
        );
    )

    boxToCell
    (
        box (0 0 0) (0 0 1);
        fieldValues
        (
            volScalarFieldValue T 600
        );
    )

Use cylinder to cell - scalar is one in circle, zero outside

::

    cylinderToCell
    (
        p1 (5 5 0);
        p2 (5 5 1);
        radius 1;
        fieldValues
        (
            volScalarFieldValue T 1
        );
    )

Example - Blob 2 (OpenFOAM 2.4.0)
---------------------------------

Copy the shockTube case

::

    $ cp -rf /home/apr207/OpenFOAM/OpenFOAM-2.4.0/tutorials/compressible/sonicFoam/laminar/shockTube/ .
    $ mv shockTube 026_blob_2
    $ cd 026_blob_2
    
Edit BCs

::

    $ nano 0/T
    
        sides
        {
            type            inletOutlet;
            inletValue      uniform 0;
            outletValue     uniform 0;
        }

    $ nano 0/p

        internalField   uniform (1 0 0);

Copy transport properties from pitzDaily tutorial (scalar transport case) and remove others

::

    $ cp $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/constant/transportProperties 
    $ rm -rf thermophysicalProperties turbulenceProperties

Set a zero diffusion

::

    $ nano constant/transportProperties
        
        DT              DT [ 0 2 -1 0 0 0 0 ] 0.0;

Change the blockMeshDict

::

    $ nano constant/polyMesh/blockMeshDict
    
        vertices
        (
            (-5 -5 -1)
            (5 -5 -1)
            (5 5 -1)
            (-5 5 -1)
            (-5 -5 1)
            (5 -5 1)
            (5 5 1)
            (-5 5 1)
        );

        blocks
        (
            hex (0 1 2 3 4 5 6 7) (100 100 1) simpleGrading (1 1 1)
        );

        edges
        (
        );

        boundary
        (
            sides
            {
                type patch;
                faces
                (
                    (1 2 6 5)
                    (0 4 7 3)
                    (3 7 6 2)
                    (0 1 5 4)
                );
            }
            empty
            {
                type empty;
                faces
                (
                    (5 6 7 4)
                    (0 3 2 1)
                );
            }
        );

Delete fvSolution, fvSchemes and controlDict (as they were for sonicFoam) and copy the ones from pitzDaily

::

    $ rm -rf controlDict  fvSchemes  fvSolution
    $ cp -rf $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/system/controlDict system/
    $ cp -rf $FOAM_TUTORIALS/basic/scalarTransportFoam/pitzDaily/system/fvS* system/

::

    $ nano system/controlDict
    
        startFrom       latestTime;

        startTime       0;

        stopAt          endTime;

        endTime         4;

        deltaT          0.01;

        writeControl    runTime;

        writeInterval   1;


::

    $ nano system/setFieldsDict

        defaultFieldValues ( volScalarFieldValue T 0.0 );

        regions          
        (
                sphereToCell
                {
                        centre (0.0 0.0 0.0);
                        radius 1;
                        fieldValues
                        (
                                volScalarFieldValue T 1
                        );
                }
        );


Create Mesh and set fields and run case

::

    $ blockMesh
    $ setFields
    $ scalarTransportFoam
    
Change the direction of the blob  - phi depends on velocity (createPhi.H), so must be deleted if velocity is changed

::

    $ rm 4/phi
    $ nano 4/U
        internalField uniform (0 1 0);
    $ nano system/controlDict
        endTime 8;
    
Change the direction of the blob  - phi depends on velocity (createPhi.H), so must be deleted if velocity is changed

::

    $ rm 8/phi
    $ nano 8/U
        internalField uniform (-1 -1 0);
    $ nano system/controlDict
        endTime 12;
        
        
