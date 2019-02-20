==============================================================
OpenFOAM Tutorials - University of Exeter (Password Protected)
==============================================================

This page contains OpenFOAM tutorials for version 2.4.0 from `Professor Gavin Tabor <http://emps.exeter.ac.uk/engineering/staff/grtabor>`_ at the University of Exeter.

.. contents::
   :local:

Tutorial 1: Basics
------------------

Download tutorial1.zip case file from VLE and save to run directory:

::

    http://vle.exeter.ac.uk/course/view.php?id=4781
    Click case files next to Tutorial 1  (introduction) - case files
    Open with Archive Manager
    Extract to /home/links/apr207/OpenFOAM/apr207-2.4.0/run

Cavity case
~~~~~~~~~~~
    
Load 2.4 and go to the tutorial1 directory and run blockMesh:

::

    of24
    run
    cd tutorial1/cavity
    blockMesh
    
Run icoFoam and foamCalc mag U to generate the magnitude of velocity at each timestep:

::

    icoFoam
    foamCalc mag U
    
Q.I.1 Using paraFoam generate plots of pressure, velocity magnitude and velocity vectors on the mid-plane of the simulation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Open the file:

::

    paraFoam
    Click Open icon
    Select cavity.OpenFOAM and press OK
    Press Apply
    Press right arrow to switch to 0.5 seconds
    Change view to 2D

Create slice:

::

    Filters > Search > Slice
    Origin: 0.05 0.05 0.005
    Normal: 0 0 1

a) Contour plot of pressure:

::

    Select color by interpolated p (dot next to p), other one is cell values
    Tools > Lock View Size Custom > 1000, 500
    Rescale to custom range -5 to 5
    Orientation Axes Visibility > Click the X
    Choose preset > Blue to Red Rainbow > Apply
    Edit color legend > Title > Pressure (m^2/s^2), Label format: %.1f > ok, Font = 9 
    Click on cavity.OpenFOAM to remove border round slice
    File > Save Screenshot > pressure
    
b) Contour plot of velocity magnitude

::

    Click slice
    Switch to interpolated magU
    Edit color legend > Title > Velocity  Magnitude (m/s), Label format: %.1f > ok, Font = 9 
    Choose preset > Blue to Red Rainbow > Apply
    Click on cavity.OpenFOAM to remove border round slice
    File > Save Screenshot > velocity
    
c) Vector plot of velocity

::

    Hide slice
    Show cavity.OpenFOAM
    Filters > Search > Cell Centers > Apply
    Filters > Search > Glyph > Apply
    Hide Cell Centers
    Glyph type > 2D Glyph
    Coloring > magU
    Scale factor: 0.005
    Turn on foam.foam > Representation > Outline
    Coloring > Solid Color > Edit > Black
    Rescale to custom range 0 to 1
    Edit color legend > Title > Velocity  Magnitude (m/s), Label format: %.1f > ok
    File > Save Screenshot > vectors

d) Save State

::

    File > Save State > Contours_and_Vectors.pvsm
    
    
QI.2 Sample dictionary to find pressure along lid and across mid-line
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

a) Run sample

::

    gedit system/sampleDict

    setFormat raw;
    
    surfaceFormat vtk;

    interpolationScheme cellPoint;

    fields
    (
        p
    );

    sets
    (    
        lineX1
        {
            type        uniform;
            axis        distance;

            start       (0.0 0.1 0.005);
            end         (0.1 0.1 0.005);
            nPoints     100;
        }
        
        lineX2
        {
            type        uniform;
            axis        distance;

            start       (0.0 0.05 0.005);
            end         (0.1 0.05 0.005);
            nPoints     100;
        }
        
    );
    
    sample

b) Plot in xmgrace

::

    xmgrace
    Data > Import > ASCII > lineX1_p.xy > OK lineX2_p.xy > OK
    Close
    Plot > Axes Properties > Xaxis > x (m) > Apply
                             Yaxis > p (Pa) > Apply
    Plot > Set Appearance > String: y = 0.1m
                            String y = 0.05m
                            
    File > Save as > plot > pressure.agr
    
    File > Print setup > png
    
    File > Print
    
QI.3 What is Re?
++++++++++++++++

Re = U*L/nu = 1*0.1/0.01 = 10

QI.4 Re-run for Re=10000
++++++++++++++++++++++++

a) Copy case to new name turbCavity

::

    cd ..
    cp -rf cavity turbCavity
    cd turbCavity
    rm -rf *.*
    
b) Change U to 1000 or nu to 0.00001, nu needs to be 1000 times less for 1000 times the Re

::

    nano constant/transportProperties

    transportModel  Newtonian;
    nu              nu [0 2 -1 0 0 0 0] 1e-05;

c) Turn on turbulence    
    
::

    nano constant/RASProperties
    
    RASModel        kEpsilon;

    turbulence      on;

    printCoeffs     on;
    
    
Run with pisoFoam

::

    nano system/controlDict

    application pisoFoam;
    
Run pisoFoam and foamCalc mag U to generate the magnitude of velocity at each timestep:

::

    pisoFoam
    foamCalc mag U    
    
How have the results changed - go back to Q.I.1 and Q.I.2
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Cylinder case
~~~~~~~~~~~~~

Go to the tutorial1 directory and run blockMesh:

::

    cd ../cylinder
    fluentMeshToFoam cylinder.msh
    
Change inlet velocity to 0.5m/s

::

    nano 0/U
    
    inlet
    {
        type            fixedValue;
        value           uniform (0.5 0 0);
    }
    
k from formula

.. math:: k = {3 \over 2} (U_{ref} T_i)^2 = {3 \over 2} (0.5 \times 0.5)^2 = 0.09375 m^2/s^2
   
::

    nano 0/k
    
    inlet
    {
        type            fixedValue;
        value           uniform 0.09375;
    }
    
Epsilon from formula

.. math:: \ell = 0.07L = 0.07 \times 5 = 0.35 m

.. math:: \varepsilon = C_{\mu}^{3/4} {{k^{3/2}} \over \ell} = 0.09^{3/4} {{0.09375^{3/2}} \over 0.35} = 0.01348 m^2/s^3
  

::
    
    nano 0/epsilon
    
    inlet
    {
        type            fixedValue;
        value           uniform 0.01348;
    }



Boundary names are in constant/polyMesh/boundary    
    
Check BCS change symmetryPlane to symmetry (symmetryPlane is for orthogonal planes, these planes are non-planar):

::

    nano 0/U 
    nano 0/p 
    nano 0/epsilon 
    nano 0/k
    nano 0/nut
    
        sides
        {
            type            symmetry;
        }

    

    
Change frontAndBackPlanes in 0/U and 0/p

::

    frontAndBackPlanes    
    {
        type            empty;
    }
    
Change schemes to bounded:

::

    nano system/fvSchemes
        
    divSchemes
    {
        default         none;
        div(phi,U)      bounded Gauss upwind;
        div(phi,k)      bounded Gauss upwind;
        div(phi,epsilon) bounded Gauss upwind;
        div(phi,R)      bounded Gauss upwind;
        div(R)          bounded Gauss linear;
        div(phi,nuTilda) bounded Gauss upwind;
        div((nuEff*dev(T(grad(U))))) Gauss linear;
    }
    
Run using simpleFoam

::

    simpleFoam

Q.I.5 Generate plots - go back to Q.I.1
+++++++++++++++++++++++++++++++++++++++

Q.I.6 Generate pressure along midline
+++++++++++++++++++++++++++++++++++++

a) CheckMesh to find domain size:

::

    checkMesh
    
    Overall domain bounding box (0 0 -0.206155) (20 5 0.206155)

a) Run sample

::

    nano system/sampleDict

    setFormat raw;

    interpolationScheme cellPoint;

    surfaceFormat vtk;
    
    fields
    (
        p
    );

    sets
    (    
        lineX1
        {
            type        uniform;
            axis        distance;

            start       (0.0 2.5 0);
            end         (5.0 2.5 0);
            nPoints     1000;
        }

        lineX2
        {
            type        uniform;
            axis        distance;

            start       (5.0 2.5 0);
            end         (20.0 2.5 0);
            nPoints     1000;
        }

    );
    
    sample


Tutorial 2: SnappyHexMesh
-------------------------

Download from VLE (ECMM148):

::

    http://vle.exeter.ac.uk/course/view.php?id=4781    
    Click case files next to Tutorial 2  (snappyhexmesh) - case files
    Open with Archive Manager
    Extract to /home/links/apr207/OpenFOAM/apr207-2.4.0/run

::

	of24
    cd shuttleSnappy
    blockMesh
    paraFoam
    File > Open > constant/triSurface/StarTrekShuttle3d.stl > Apply
    Exit paraview

Run snappyHexMesh

::

    snappyHexMesh
    checkMesh

It will say there is an error

Max skewness = 4.694484, 2 highly skew faces detected which may impair the quality of the results

::

    foamToVTK -faceSet skewFaces -time 3

Open Paraview and view the skewed faces:

::

    File > Open VTK/skewFaces
    Color them red

Correct the problem with skew faces:

Max refinement surfaces from 4 to 5

::

    refinementSurfaces
    {
        starTrekShuttle
        {
            // Surface-wise min and max refinement level
            level (2 5);
        }
    }

::

    rm -rf 1 2 3

::

    snappyHexMesh
    checkMesh

::

    paraFoam
    Move to time = 3
    Click on Clip
    Apply
    Surface with edges
    Filters > Alphabetical > Extract Cells by Region > Apply
    Filters > Alphabetical > Extract Block

    Show only ExtractBlock1
    Show internalMesh and zoom in to view

Mesh is ok    
Delete polyMesh folder in shuttleSmall
Copy folder 3/polyMesh to shuttleSmall/constant/polyMesh

From the over directory:

::

    rm -rf shuttleSmall/constant/polyMesh
    cp -rf shuttleSnappy/3/polyMesh shuttleSmall/constant

::

    checkMesh
    transformPoints -scale "(0.01 0.01 0.01)"
    checkMesh

Change constant/polyMesh/boundary starTrekShuttle to starTrekShuttle_AutoCAD

::

    gedit 0/* constant/polyMesh/bundary
    rm -rf shuttleSmall/0
    cp -rf shuttleSmall/0.orig shuttleSmall/0

::

    starTrekShuttle_AutoCAD


::

    cp -rf 0.orig/* 0

::

    potentialFoam

--> FOAM FATAL IO ERROR: 
keyword laplacian(1,Phi) is undefined in dictionary "/home/apr207/OpenFOAM/apr207-2.4.0/run/tutorial2/shuttleSmall/system/fvSchemes.laplacianSchemes"

::

    gedit system/fvSchemes
    laplacian(1,Phi) Gauss linear corrected;

::

    potentialFoam

--> FOAM FATAL IO ERROR: 
keyword Phi is undefined in dictionary "/home/apr207/OpenFOAM/apr207-2.4.0/run/tutorial2/shuttleSmall/system/fvSolution.solvers"

Copy the solver for p to Phi

::

    Phi
    {
    	solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0;
        minIter         0;
    	maxIter         1000;
    }

::

    potentialFoam

--> FOAM FATAL ERROR : flux requested but Phi not specified in the fluxRequired sub-dictionary of fvSchemes. 

::

    fluxRequired
    {
        default         no;
        p;
        Phi;
    }


Run simpleFoam - need to remove calcForces so that the post processing works.

::

    rm -rf calcForces
    simpleFoam >& log &
    tail -f log

::

    foamLog log


Can plot .xy in xmgrace for UxFinalRes_0, UyFinalRes_0, UzFinalRes_0, epsilonFinalRes_0, kFinalRes_0, pFinalRes_1

Can load in Excel foamCalc forces, use spaces and use ((())) as delimiters.

Servers:

::

    scp -r shuttleSnappy apr207@emps-ugeng1:
    ssh userID@emps-ugeng1
    nohup nice -19 timeout 1h simpleFoam > log &

To find out the load on the servers:

::

    uptime

example result:
13:51:55 up 77 days,  3:20,  3 users,  load average: 0.07, 0.06, 0.05
The first number after "load average" is the most relevant. An average of 1 represents 1 core running at 100\%. So if it reaches 28 the server is running at full capacity.

Blue or purple room:

::

    nohup nice -19 simpleFoam > log &
    logout



Tutorial 3: Backward Facing Step
--------------------------------

Preparation
~~~~~~~~~~~

Changed inlet height to 2H from H - changes total height of domain to 3:

Origin is located at base of step

::

    of24
    run
    cd bfStep
    nano constant/polyMesh/blockMeshDict

    vertices        
    (
        (0 0 0)
        (10 0 0)
        (10 1 0)
        (0 1 0)
        (0 0 0.1)
        (10 0 0.1)
        (10 1 0.1)
        (0 1 0.1)

        (10 3 0)
        (0 3 0)
        (10 3 0.1)
        (0 3 0.1)

        (-2 1 0)
        (-2 3 0)
        (-2 1 0.1)
        (-2 3 0.1)
    );

Blocks 2 and 3 have doubled in height, so double number of cells    
    
::
    
    blocks          
    (
        hex (0 1 2 3 4 5 6 7) (40 10 1) simpleGrading (1 1 1)
        hex (3 2 8 9 7 6 10 11) (40 20 1) simpleGrading (1 1 1)
        hex (12 3 9 13 14 7 11 15) (8 20 1) simpleGrading (1 1 1)
    );
    
Run blockMesh

::

    blockMesh
    
Run checkMesh

::

    checkMesh

Modified fluid viscosity to achievve Re = 44,000 = U 2H / nu, nu = 5 x 2 x 0.1 / 44000 = 2.27e-5    

::

    nano constant/transportProperties
    
    transportModel  Newtonian;

    nu              nu [0 2 -1 0 0 0 0] 2.27e-05;

Modified k and epsilon according to an inlet length scale of 2H (0.2) and an inlet turbulent intensity of 5% (0.05). k = 0.09375 m^2/s^2. epsilon = 0.337. 

::

    nano 0/k
    
    internalField   uniform 0.09375;

    boundaryField
    {
        top
        {
            type            kqRWallFunction;
            value           uniform 0.09375;
        }
        step
        {
            type            kqRWallFunction;
            value           uniform 0.09375;
        }
        base
        {
            type            kqRWallFunction;
            value           uniform 0.09375;
        }
        inlet
        {
            type            fixedValue;
            value           uniform 0.09375;
        }
        outlet
        {
            type            zeroGradient;
        }
        frontAndBack
        {
            type            empty;
        }
    }

::

    nano 0/epsilon
    
    internalField   uniform 0.337;

    boundaryField
    {
        top
        {
            type            epsilonWallFunction;
            Cmu             0.09;
            kappa           0.41;
            E               9.8;
            value           uniform 0.337;
        }
        step
        {
            type            epsilonWallFunction;
            Cmu             0.09;
            kappa           0.41;
            E               9.8;
            value           uniform 0.337;
        }
        base
        {
            type            epsilonWallFunction;
            Cmu             0.09;
            kappa           0.41;
            E               9.8;
            value           uniform 0.337;
        }
        inlet
        {
            type            fixedValue;
            value           uniform 0.337;
        }
        outlet
        {
            type            zeroGradient;
        }
        frontAndBack
        {
            type            empty;
        }
    }
    
    
Starting Point - kEpsilon model - four mistakes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Change a mistake (although it runs) in system/fvSolution - remove PBiCG from k entry

::

    nano system/fvSolution
    
    k
    {
	solver		 PBiCG;
        preconditioner   DILU;
        tolerance        1e-05;
        relTol           0.1;
    };

Another mistake: change turbulenceModel kEpsilon; to simulationType RASModel;

::

    nano constant/turbulenceProperties
    simulationType RASModel;

    
Another mistake: Remove UyFinal from system  folder    
    
::

    rm -rf system/UyFinal
    
Another mistake in control Dict, turbFoam isn't an application    
    
::

    nano system/controlDict
    application     simpleFoam;

Better starting point
~~~~~~~~~~~~~~~~~~~~~

Copy bfStep to kEpsilon and run simpleFoam on kEpsilon model

::

    cp -rf bfStep kEpsilon
    cd kEpsilon
    simpleFoam > log 2>&1 &

View the log file using less

::

    less log

Go to the end of the log file:

::

    Shift G

Show that Initial residual for p didn't converge until 160 iterations

Quit the log file

::

    q

kEpsilon RNG
~~~~~~~~~~~~

::

    cd ..
    cp -rf bfStep RNGkEpsilon
    cd RNGkEpsilon
    

Change turbulence model to RNGkEpsilon

::

    nano constant/RASProperties
    turbulenceModel RNGkEpsilon;
    
    
    
Run simpleFoam

::

    simpleFoam > log 2>&1 &
    
Use less

::

    less log

Go to end:

::

    shift g
    
Go to start:

::

    g


    
kOmega
~~~~~~

::

    cd ..
    cp -rf bfStep kOmega
    cd kOmega

Change turbulence model to kOmega

::

    nano constant/RASProperties
    RASModel        kOmegaSST;  
    
Create omega file    

::

    cp -rf 0/epsilon 0/omega
    nano 0/omega

Omega = epsilon/k = 0.337 / 0.09375 = 3.594   

Units = m^2/s^2 / m^2/s^3 = s^-1
    
::

    dimensions      [0 0 -1 0 0 0 0];

    internalField   uniform 3.594;

    boundaryField
    {
        top
        {
            type            omegaWallFunction;
            Cmu             0.09;
            kappa           0.41;
            E               9.8;
            value           uniform 3.594;
        }
        step
        {
            type            omegaWallFunction;
            Cmu             0.09;
            kappa           0.41;
            E               9.8;
            value           uniform 3.594;
        }
        base
        {
            type            omegaWallFunction;
            Cmu             0.09;
            kappa           0.41;
            E               9.8;
            value           uniform 3.594;
        }
        inlet
        {
            type            fixedValue;
            value           uniform 3.594;
        }
        outlet
        {
            type            zeroGradient;
        }
        frontAndBack
        {
            type            empty;
        }
    }

open fvSchemes to copy divSchemes and laplacianSchemes for omega:

::

    nano system/fvSchemes
    
    divSchemes
    {
        div(phi,omega) Gauss upwind;
    }
    
    laplacianSchemes
    {
        laplacian(DomegaEff,omega) Gauss linear corrected;
    }

open system/fvSolution

::

    nano system/fvSolution
    
    omega
    {
        solver           PBiCG;
        preconditioner   DILU;
        tolerance        1e-05;
        relTol           0.1;
    };
    
    relaxationFactors
    {
        omega           0.7;
    }

Run simpleFoam

::

    simpleFoam > log 2>&1 &

     
LRR Reynolds Stress - start from stratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    cd ..
    cp -rf bfStep LRR
    cd LRR
    
Run R utility

::

    R
    
Change turbulence model

::

    nano constant/RASProperties
    
    RASModel LRR;
    
Run simpleFoam

::

    simpleFoam > log 2>&1 &   
    

LRR Reynolds Stress - start from kEpsilon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    cd ..
    cp -rf kEpsilon LRR2
    cd LRR2
    
Run R utility

::

    R
    
Change turbulence model

::

    nano constant/RASProperties
    
    RASModel LRR;

Start from latest time:

::

    nano system/controlDict

    startFrom       latestTime;

    
Run simpleFoam

::

    simpleFoam > log 2>&1 &     

Run wallGradU utility    
    
::

    wallGradU
    
    
Plot streamlines using SurfaceLIC

::

    paraFoam
    
Tools > Manage Plugins > SurfaceLIC > AutoLoad    

File > Exit

::

    paraFoam
    
Change to last time    
    
Change to U

Change to SurfaceLIC

Plot wallGradU_X
~~~~~~~~~~~~~~~~

Select header and select all elements > Apply

Filter > Extract Block > select base > Apply

Show header again

Switch back to ExtractBlock1

Filters > Plot Data

Select wallGradU_X or wallGradU(0)

Check with a visual, switch to header and select streamline screen

Filters >  Plot over line

Enter:

0.43 0 0
0.43 0 0

Save state and load in kEpsilon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File > Save State > state

Close paraview

::

    cd ..

    cp -rf LRR2/state.pvsm kEpsilon/state.pvsm 

    cd kEpsilon

    wallGradU

paraview

File > Load State

Select state

Change .OpenFOAM file location to kEpsilon

Refresh (if needed)

Plot velocity profile around 0.4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Filters > Plot Over Line

0.4 0.0 0.005
0.4 0.3 0.005

Deselect everything except U_X

Find point were there is no reversed flow (negative velocity)

Move to 0.425 (shows -ve velocity)

Move to 0.426 (shows no -ve velocity)


Tutorial 4: Pointwise
---------------------

Preparation
~~~~~~~~~~~

Download NACA6412.igs geometry from VLE   

Open pointwise

::

    pointwise
    
File > Import > Database

Select NACA6412.igs from browser > Open > OK

Defaults > Dimension = 75

List > Database > Select all four database entries > Click "Connectors on Database Entries" in toolbar

List > Connectors > Select all connectors > Points on

Uncheck the database mask (to put database in grey

Zoom in to trailing edges and select both

Grid > Dimension > Average delta s = 0.0005 > Dimension > OK

Click "All Masks On/Off"

Check the Spacing Constaints Mask

Select spacing constraints on upper and lower surfaces of the aerofoil at the trailing edge

Click spacing on the toolbar and enter 0.0005

Select spacing constaints on the upper and lower surace of the aerofoil at the leading edge

Spacing = 0.005

Click Connectors mask and select all connectors

Create > Extrude > Normal > Done

Attributes tab > Step Size frame > Initial delta s = 0.0001

               > Orientation frame > Flip

Run tab > 91 steps > Run > OK

Click Save on toolbar > aerofoil_after_pointwise_tutorial.pw

Download OpenFOAM case from ele
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download Tutorial 4

Mesh generation - Gavin tutorial
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

---

**Create a 3D block**

Open aerofoil.pw

CAE > Select Solver > OpenFOAM

Click "All Masks On/Off" to mask all entities and unmask domains

Select dom-1

Create > Extrude > Translate

Steps = 1

Direction = 0.0 0.0 1.0

Distance = 0.25

Run 

OK

---

CAE > Set boundary conditions

New > Name = farfield > CAE Type = patch

View > Set view type > Perspective

Control + Right mouse = rotate view

Select the inlet > click the tick

Click off into the blue

New > Name = frontAndBackPlanes > CAE Type = empty

Select the front and back planes > click the tick

New > Name = airfoil > CAE Type = wall

Click on boundary 0 > click Add to Selection

Close

---

File > Save

Click "All Masks On/Off"

Select blocks blk-1

File > Export > CAE > Navigate to tutorial4/NACA6412/constant/ Create polyMesh folder > Select it > Choose

Navigate to the case NACA6412 directory and run paraFoam to check mesh

---


Renumber mesh to enchance diagonal dominance:

::

    renumberMesh -overwrite

It will complain that the type of front and back planes are not defined:

::
    
    frontAndBackPlanes
    {
        type empty;
    }   

Add this to controlDict:    
    
::

    writeInterval   10;


    functions
    {
        calcForces
        {
            type                forces;
            functionObjectLibs  ( "libforces.so" );
            enabled             true;
            outputControl       outputTime;
            pName               p;
            UName               U;
            rhoName             rhoInf;
            rhoInf              1.2;
            CofR                (0.25 0 0);
            patches
            (
                airfoil
            );
        }
    }   
    
    
Run simpleFoam    
    
::

    rm -rf forces
    simpleFoam > log 2>&1 & 

::

    foamLog log
    
    
Mesh for Pipe Flow - Use Pointwise Tutorial
+++++++++++++++++++++++++++++++++++++++++++

http://www.pointwise.com/DIY/Butterfly-Topology.shtml

Create a circle with diameter 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Script > Begin journaling > script.glf

Create > Draw Curves > Circle > Entity type = database > Segment type = Circle 

> XYZ = 0 0 0 > Press enter

> XYZ = 1 0 0 > Press enter

2 Points and angle > Angle = 360.0 > Press enter > OK

Split the circle into four
^^^^^^^^^^^^^^^^^^^^^^^^^^

Select the circle

Edit > Split > Persent of Length = 50 > Press enter > OK (Repeat creating four segments)

Defaults > Dimension = 20 > Press enter

Create domain
^^^^^^^^^^^^^

Select all databases > Connectors on Database Entities

Select all connectors > Assemble Domains

Examine domain
^^^^^^^^^^^^^^

Examine > Maximum included angle (shows 175 degree angle)

Select domain > Press delete key

Create curves
^^^^^^^^^^^^^

Create > 2 Point curves

Click top and bottom points

Click left and right points

Select horizontal connector > Edit > Split > Percent of Length = 50 > Press enter > OK (Repeat for vertical connector and radii)

Create > 2 Point curves (Create central box)

Remove interior lines within box

Select all connectors > Dimension = 20

Select quadrants > Assemble domains

Examine domain
^^^^^^^^^^^^^^

Examine > Maximum included angle (shows 135 degree angle)

Select domain > Press delete key

Solver
^^^^^^

Select domains > Grid > solve

Edge Attributes > Select all interior edges (twice) by clicking a box select (left click) over the connectors

Boundary Conditions > Floating > Apply

Solve > Run > Iterations = 100 > OK

Examine domain
^^^^^^^^^^^^^^

Examine > Maximum included angle (shows 122 degree angle)

Extrude domain
^^^^^^^^^^^^^^

CAE > Select Solver > OpenFOAM

Select all domains

Create > Extrude > Translate > Done

Steps = 100

Direction = 0.0 0.0 1.0

Distance = 10

(Creates a domain with 200,000 cells)

Boundary conditions
^^^^^^^^^^^^^^^^^^^

CAE > Set boundary conditions

New > Name = inlet > CAE Type = patch

Select the inlet > click the tick

Click off into the blue

New > Name = outlet > CAE Type = patch

Select the outlet > click the tick

New > Name = wall > CAE Type = wall

Click on boundary 0 > click Add to Selection

Close

Save
^^^^

File > Save

Script > End journaling

Export
^^^^^^

Click “All Masks On/Off”

Select blocks blk-1

File > Export > CAE > Navigate to case_directory/constant/ Create polyMesh folder > Select it > Choose

Navigate to the case directory and run checkMesh and paraFoam to check mesh
  
    
Tutorial 5 - Free Surface Flow
------------------------------

Download case to andrew-2.4.0/run/

Test swak4Foam is installed:

::

    two
    funkySetFields

Change to weir case directory and run blockMesh    
    
::

    cd weir
    blockMesh

interFoam factors out the hydrostatic pressure: :math:`p^* = p - \rho g y`

Tutorial 6 - Heat Transfer
--------------------------

Download case to apr-2.4.0/run/

If you use the solver, then:

/home/andrew/OpenFOAM/OpenFOAM-2.4.0/src/finiteVolume/lnInclude/cyclicAMIFvPatch.H:39:35: fatal error: cyclicAMILduInterface.H: No such file or directory
compilation terminated.
  
    
User directory
~~~~~~~~~~~~~~
    
What are applications?
++++++++++++++++++++++

Executables - has a main() function.

For applications (executables), e.g. solvers, the user location is:

::

    $WM_PROJECT_USER_DIR/applications/solvers

What are libraries?
+++++++++++++++++++

Run-time selectable classes - doesn't have a main() function.

For library source (e.g. viscosity models):

::

    $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
    
    
Compilation
~~~~~~~~~~~

What is wclean?
+++++++++++++++

Cleans directory of dependency files (run before wmake)

::

    wclean

What is wmake?
++++++++++++++

Compilation script based on make.


What does wmake do?
+++++++++++++++++++

It compiles OpenFOAM code in the source code directory. 

How is wmake run?
+++++++++++++++++

In the source code directory (which can be placed anywhere), although usually in the **user directory**. Contains at least:

- solver.C file
- Make directory containing files and options files

Run using:

::

    wmake
    
What is rehash?
+++++++++++++++

Informs compiler that something has been updated (run after wmake). Makes the executable available. NOT USED ANYMORE!

::

    rehash   
    

Files and options
~~~~~~~~~~~~~~~~~

What is the Make/files file?
++++++++++++++++++++++++++++

Specifies:

- Which files need to be complied (cassonFoam.C)
- Where the resulting executable will be placed and what it's name is (EXE = $(FOAM_USER_APPBIN)/cassonFoam)

::

    cassonFoam.C

    EXE = $(FOAM_USER_APPBIN)/cassonFoam

- Shouldn't write to FOAM_APPBIN (installation direcctory), instead write to FOAM_USER_APPBIN (user directory)        
- File name need not be the same at the executable name
- File name need not be unique
- Exectuable name has to be unique

Download from ELE
~~~~~~~~~~~~~~~~~

Extract to andrew-2.4.0/run/Project/

Copy boussinesqFoam to andrew-2.4.0/applications/solvers

::

    mv boussinesqFoam/ $WM_PROJECT_USER_DIR/applications/solvers/

Edit createFields.H
~~~~~~~~~~~~~~~~~~~    
    
- Read additional properties
- Add a field for theta

::

    Info<< "Reading field theta\n" << endl;
    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    dimensionedScalar kappa
    (
        transportProperties.lookup("kappa")
    );
    
    dimensionedScalar rho0
    (
        transportProperties.lookup("rho0")
    );

    dimensionedScalar Cv
    (
        transportProperties.lookup("Cv")
    );
    
    dimensionedScalar theta0
    (
        transportProperties.lookup("theta0")
    );    

    dimensionedScalar beta
    (
        transportProperties.lookup("beta")
    ); 

    dimensionedScalar hCoeff = kappa/(rho0*Cv);
    
    dimensionedVector g
    (
        transportProperties.lookup("g")
    );

    
Edit icoFoam.C
~~~~~~~~~~~~~~

- Add source term to momentum equation
- Add temperature equation

::

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
        + fvm::div(phi, U)
        - fvm::laplacian(nu, U)
        + beta*g*(theta0-theta)
    );

    solve(UEqn == -fvc::grad(p));

Where to add the temperature equation?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Coupling from momentum to temperature equation due to convection
- No coupling from temperature to momentum equation
- Temperature equation can be added at the end of each timestep

::

    fvScalarMatrix tempEqn
    (
        fvm::ddt(theta)
        + fvm::div(phi, theta)
        - fvm::laplacian(hCoeff, theta)
    );

    tempEqn.solve();        
    
    runTime.write();

Compile
~~~~~~~

::

    wclean
    wmake    
    
/home/andrew/OpenFOAM/OpenFOAM-2.4.0/src/finiteVolume/lnInclude/cyclicAMIFvPatch.H:39:35: fatal error: cyclicAMILduInterface.H: No such file or directory
compilation terminated.

Add meshTools:

::

    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude

    EXE_LIBS = \
        -lfiniteVolume \
        -lmeshTools    
    

icoFoam.C:74:19: error: ‘ddtPhiCorr’ is not a member of ‘Foam::fvc’
                 + fvc::ddtPhiCorr(rAU, U, phi);

PISO Loop has completely changed in 2.4.0

::
                 
    // --- PISO loop

    for (int corr=0; corr<nCorr; corr++)
    {
        volScalarField rAU(1.0/UEqn.A());

        volVectorField HbyA("HbyA", U);
        HbyA = rAU*UEqn.H();
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            (fvc::interpolate(HbyA) & mesh.Sf())
            + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
        );

        adjustPhi(phiHbyA, U, p);

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi = phiHbyA - pEqn.flux();
            }
        }

        #include "continuityErrs.H"

        U = HbyA - rAU*fvc::grad(p);
        U.correctBoundaryConditions();
    }    

Create the mesh
~~~~~~~~~~~~~~~

::

    /OpenFOAM/andrew-2.4.0/run/Project/tutorial6/bendHeat
    blockMesh
    
Create the dictionaries
~~~~~~~~~~~~~~~~~~~~~~~

::

    cp -rf 0/U 0/theta
    nano theta
        FoamFile
        {
            version     2.0;
            format      ascii;
            class       volScalarField;
            object      theta;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
        dimensions       [0 0 0 1 0 0 0];

        internalField    uniform 300;

        boundaryField
        {
            inlet
            {
                type                fixedValue;
                value               uniform 300;
            }

            outlet
            {
                type                zeroGradient;
            }

            walls
            {
                type                fixedValue;
                value               uniform 350;
            }

            frontAndBack
            {
                type                empty;
            }
        }

Check transportProperties

::

    nano constant/transportProperties
        
        // Water values
        nu              nu [0 2 -1 0 0 0 0] 1.0e-6;

        Cv              Cv [0 2 -2 -1 0 0 0] 4182;

        rho0            rho0 [1 -3 0 0 0 0 0] 998;

        kappa           kappa [1 1 -3 -1 0 0 0] 0.6;

        beta            beta [0 0 0 -1 0 0 0] 0.00021;

        theta0          theta0 [0 0 0 1 0 0 0] 300;

        g               g [0 1 -2 0 0 0 0] (0.0 -9.81 0.0);    

Check fvSchemes
        
::

    nano system/fvSchemes
    
    ddtSchemes
    {
        default Euler;
    }

    gradSchemes
    {
        default         Gauss linear;
        grad(p)         Gauss linear;
    //    grad(theta)         Gauss linear;
    }

    divSchemes
    {
        default         none;
        div(phi,U)      Gauss linear;
        div(phi,theta)  Gauss linear;
    }

    laplacianSchemes
    {
        default         Gauss linear corrected;
        laplacian(nu,U) Gauss linear corrected;
        laplacian((1|A(U)),p) Gauss linear corrected;
    //    laplacian(kappa/(rho0*Cv),theta) Gauss linear corrected;
    }

    interpolationSchemes
    {
        default         linear;
        interpolate(HbyA) linear;
    }

    snGradSchemes
    {
        default         corrected;
    }

    fluxRequired
    {
        default         no;
        p;
    }
   
Check fvSolution
   
::

    nano system/fvSolution
    
    solvers
    {
        p
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-06;
            relTol          0;
        }

        U
        {
            solver          PBiCG;
            preconditioner  DILU;
            tolerance       1e-05;
            relTol          0;
        }

        theta
        {
            solver          PBiCG;
            preconditioner  DILU;
            tolerance       1e-05;
            relTol          0;
        }
    }

    PISO
    {
        nCorrectors     2;
        nNonOrthogonalCorrectors 0;
        pRefCell        0;
        pRefValue       0;
    }

Run the solver
~~~~~~~~~~~~~~

::

    boussinesqFoam
    
Evaluate Re, Pr and Nu
~~~~~~~~~~~~~~~~~~~~~~

.. math::

    Re = {{U D} \over \nu} = {{0.01 \times 0.01} \over {1.0 \times 10^{-6} }} = 100

.. math::

    Pr = {{c_p \mu} \over k} = {{c_p \nu \rho} \over k} = {{4182 \times 1.0 \times 10^{-6} \times 998} \over 0.6} = 7

.. math::

    Nu = {{h L} \over k} = {{{k \over { \rho c_v}} L} \over k} = {L \over { \rho c_v}} = {0.01 \over {998 \times 4182}} = 2 \times 10^{-9} 
    
Plot Contours of Temp and Velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Shows that simulation is transient and has not reached steady state.

Plotting outlet temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Shows that it is changing

Create Average Temperature Field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In icoFoam.C

::

    double n = 0;
    while (runTime.loop())
    {
        .
        .
        .
        tempEqn.solve(); 
        thetaAv = (n / (n+1.0))*thetaAv+(1.0/(n+1.0))*theta;
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  Timesteps = " << n << " -"
            << nl << endl;
        n = n+1.0;
    }

In createFields.H    
    
::
    
    Info<< "Reading field theta average\n" << endl;
    volScalarField thetaAv
    (
        IOobject
        (
            "thetaAv",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("thetaAv", dimensionSet(0,0,0,1,0,0,0), scalar(1.0))
    );    
    
    
Project 2 - Casson Model
------------------------


Mesh for parallel plates = modify bfStep
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hydraulic Diameter
++++++++++++++++++

H = 1, so D_h = 2H = 2m


Reynolds Number for laminar flow, say 500
+++++++++++++++++++++++++++++++++++++++++

Re = U D_h / nu

U = 1

nu = (U D_h)/Re = (1 * 2)/500 = 0.004


Development length
++++++++++++++++++

Schlichting:

L = 0.16 (H/2) * (U * (H/2) / nu) = 10

Choose 20


Copy across case from offSetCylinder
++++++++++++++++++++++++++++++++++++

- Change to case directory

::

    cd $WM_PROJECT_USER_DIR/run/Project/

- Copy files

::

    cp -rf $FOAM_TUTORIALS/incompressible/nonNewtonianIcoFoam/offsetCylinder/ 002_plates_nonNewtonianIcoFoam


Mesh for Plates - Use bfStep
++++++++++++++++++++++++++++

Can use this tool to check for a cell-to-cell expansion ratio of 1.2 - must specify overall expansion ratio

http://openfoamwiki.net/index.php/Scripts/blockMesh_grading_calculation

hex blocks: 

- x1 direction is vertex 0 to 1
- x2 direction is vertex 1 to 2
- x3 plane defined by vertices 0 1 2 3

Right hand rule for defining 

::

    convertToMeters 1;

    // The vertices can be defined using variables Lx, Ly, Lz. It saves time to modify the size of the domain

    Lx 20;
    Ly 1;
    Lz 0.1;

    vertices
    (
        (0 0 0)
        ($Lx 0 0)       // First vertex in x1 direction
        ($Lx $Ly 0)     // Second vertex in x2 direction
        (0 $Ly 0)       // Third vertex completes plane
        (0 0 $Lz)
        ($Lx 0 $Lz)
        ($Lx $Ly $Lz)
        (0 $Ly $Lz)
    );

    // Mesh definition (homogeneous grid with a single cell in the y and z axis since the simulation is 1D)

    blocks
    (
        hex (0 1 2 3 4 5 6 7) (200 20 1) 
        simpleGrading 
        (
            1	// x-direction expansion
            (
                    (0.2 0.3 2.5)             // 20% y-dir, 30% cells (6 cells), expansion = 2.5
                    (0.6 0.4 1)             // 60% y-dir, 40% cells (8 cells), expansion = 1
                    (0.2 0.3 0.4)          // 20% y-dir, 30% cells (6 cells), expansion = 0.4 (1/2.5)
            )
            1       // z-direction expansion
        )
    );

    boundary
    (
        inlet
        {
            type patch;
            faces
            (
                            (0 4 7 3)
            );
        }
        outlet
        {
            type patch;
            faces
            (
                            (2 6 5 1)
            );
        }

            // Faces orthogonal to y and z axis are defined as «empty» to specify that the simulation is 1D

        walls
        {
            type wall;
            faces
            (
                (3 7 6 2)
                (0 1 5 4)
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

Boundary Conditions
+++++++++++++++++++
    
Change names of boundary conditions and remove irrelevant ones (left is inlet, right is outlet, cylinder is walls, defaultFaces is frontAndBack)
    
::
    
    nano 0/p

        dimensions      [0 2 -2 0 0 0 0];

        internalField   uniform 0;

        boundaryField
        {
            inlet
            {
                type            zeroGradient;
            }

            outlet
            {
                type            fixedValue;
                value           uniform 0;
            }

            frontAndBack
            {
                type            empty;
            }
        }

::

    nano 0/U
    
        dimensions      [0 1 -1 0 0 0 0];

        internalField   uniform (0 0 0);

        boundaryField
        {
            inlet
            {
                type            fixedValue;
                value           uniform (1 0 0);
            }

            outlet
            {
                type            zeroGradient;
            }

            walls
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }

            frontAndBack
            {
                type            empty;
            }
        }

Kinematic viscosity
+++++++++++++++++++

::

    nano constant/transportProperties
    
        transportModel  Newtonian;

        nu              nu [ 0 2 -1 0 0 0 0 ] 0.004;


Can increase timestep by factor of 10
+++++++++++++++++++++++++++++++++++++

As Courant number order of 0.03, can increase timestep:

::

    deltaT          0.025;

Also need to increase time to 20 seconds to get profile similar to analytical

::

    endTime         20;
    
Takes around 3 seconds (4,000 cells)


Mesh for cylinder = modify example on openfoamwiki
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

https://openfoamwiki.net/index.php/Main_ContribExamples/AxiSymmetric

::

    convertToMeters 1;

    Lz 20;
    Lx 0.999;
    Ly 0.044;
    Ly2 -0.044;
    
    vertices
    (
        (0 0 0)
        ($Lx $Ly 0)
        ($Lx $Ly $Lz) 
        (0 0 $Lz)
        ($Lx $Ly2 0) 
        ($Lx $Ly2 $Lz)
    );
    
    blocks
    (
    hex (0 4 1 0 3 5 2 3) (20 1 200) simpleGrading (0.03125 1 1)
    );
    
    edges
    (
    );
    
    boundary
    (
        front
        { 
            type wedge;
            faces  
            (
                (0 1 2 3)
            );
        }
        back
        { 
            type wedge;
            faces  
            (
                (0 3 5 4)
            );
        }
        pipeWall
        { 
            type wall;
            faces  
            (
                (1 4 5 2)
            );
        }
        inlet
        { 
            type patch;
            faces  
            (
                (0 4 1 0)
            );
        }
        outlet
        { 
            type patch;
            faces  
            (
                (3 5 2 3)
            );
        }
        axis
        { 
            type empty;
            faces  
            (
                (0 3 3 0)
            );
        }
    );

Boundary Conditions
+++++++++++++++++++
    
Add front, back and axis.

Change velocity direction
    
::
    
    nano 0/p

        dimensions      [0 2 -2 0 0 0 0];

        internalField   uniform 0;

        boundaryField
        {
            inlet
            {
                type            zeroGradient;
            }

            outlet
            {
                type            fixedValue;
                value           uniform 0;
            }

            frontAndBack
            {
                type            empty;
            }
            front
            {
                type            wedge;
            }

            back
            {
                type            wedge;
            }

            axis
            {
                type            empty;
            }
        }

::

    nano 0/U
    
        dimensions      [0 1 -1 0 0 0 0];

        internalField   uniform (0 0 0);

        boundaryField
        {
            inlet
            {
                type            fixedValue;
                value           uniform (0 0 1);
            }

            outlet
            {
                type            zeroGradient;
            }

            walls
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }
            front
            {
                type            wedge;
            }

            back
            {
                type            wedge;
            }

            axis
            {
                type            empty;
            }
        }
    
    
Hardcoding in the icoFoam Solver - cassonFoam
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Copy icoFoam solver to apr207-2.4.0/applications/solvers/cassonFoam

::

    cd $WM_PROJECT_USER_DIR/applications/solvers
    cp -rf $FOAM_APP/solvers/incompressible/icoFoam cassonFoam
    mv icoFoam.C cassonFoam.C
    nano Make/files
    
        cassonFoam.C

        EXE = $(FOAM_USER_APPBIN)/cassonFoam
    wclean
    wmake

Edit createFields.H

::

    // ADDED

    Info  << "Reading field nu" << nl << endl;
    volScalarField  nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("nu", dimensionSet(0,2,-1,0,0,0,0), scalar(1.0))
    );

    dimensionedScalar eta
    (
        transportProperties.lookup("eta")
    );

    dimensionedScalar tau_y
    (
     	transportProperties.lookup("tau_y")
    );

    dimensionedScalar rho
    (
     	transportProperties.lookup("rho")
    );
    
    // ADDED
    
Edit cassonFoam.C - after the PISO loop has finished:    
    
::

    // --- PISO loop

    for (int corr=0; corr<nCorr; corr++)
    {
    .
    .
    .
    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    }
    
    volScalarField  J2 = 0.5* magSqr(symm(fvc::grad(U))) + dimensionedScalar("tiny",dimensionSet (0,0,-2,0,0,0,0),0.0001);   // ADDED
    nu = sqr(pow(sqr(eta)*J2 ,0.25) + pow(tau_y /2 ,0.5))/( rho*sqrt(J2));                                                   // ADDED

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
                
transportProperties    
    
::
    
    eta        eta [ 1 -1 -1 0 0 0 0 ] 4.86;
    tau_y     tau_y [ 1 -1 -2 0 0 0 0 ] 14.38;
    rho        rho [1  -3 0 0 0 0 0]  1200;   

Run cassonIcoFoam    
    
::

    cassonFoam
    foamCalc mag U
    sample
    
    
Runtime selectable viscosity model - CassonGavin
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy BirdCarreau into directory:

::

    $ cd $WM_PROJECT_DIR
    $ cp -r --parents src/transportModels/incompressible/viscosityModels/BirdCarreau/ $WM_PROJECT_USER_DIR/
    $ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
    $ mv BirdCarreau CassonGavin
    $ cd CassonGavin
    $ mv BirdCarreau.C CassonGavin.C
    $ mv BirdCarreau.H CassonGavin.H
    $ sed -i s/BirdCarreau/CassonGavin/g CassonGavin.C
    $ sed -i s/BirdCarreau/CassonGavin/g CassonGavin.H
    $ sed -i s/BirdCarreau/CassonGavin/g Make/files
    

Add formula for CassonGavin plot - private member variables denoted by underscore e.g. $U_$

::

    $ nano CassonGavin.C
    
    #include "CassonGavin.H"
    #include "addToRunTimeSelectionTable.H"
    #include "surfaceFields.H"
    #include "fvcGrad.H"                        //ADDED
    // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    namespace Foam
    {
    namespace viscosityModels
    {
        defineTypeNameAndDebug(CassonGavin, 0);
        addToRunTimeSelectionTable
        (
            viscosityModel,
            CassonGavin,
            dictionary
        );
    }
    }
    
        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        Foam::tmp<Foam::volScalarField>
        Foam::viscosityModels::CassonGavin::calcNu() const
        {
        
            volScalarField  J2 = 0.5* magSqr(symm(fvc::grad(U_))) + dimensionedScalar("tiny",dimensionSet (0,0,-2,0,0,0,0),0.0001);   // ADDED
            return (sqr(pow(sqr(eta_)*J2 ,0.25) + pow(tau_y_ /2 ,0.5))/( rho_*sqrt(J2)));                                           // ADDED

        }

Add new coefficients:

::

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::viscosityModels::CassonGavin::CassonGavin
    (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
    )
    :
            viscosityModel(name, viscosityProperties, U, phi),
            CassonGavinCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
            eta_(CassonGavinCoeffs_.lookup("eta")),         // ADDED
            tau_y_(CassonGavinCoeffs_.lookup("tau_y")),     // ADDED
            rho_(CassonGavinCoeffs_.lookup("rho")),         // ADDED
            nu_
            (
                IOobject
                (
                    "nu",
                    U_.time().timeName(),
                    U_.db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                calcNu()
            )
    {}

Add read coefficients:

::

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    bool Foam::viscosityModels::CassonGavin::read
    (
            const dictionary& viscosityProperties
    )
    {
            viscosityModel::read(viscosityProperties);

            CassonGavinCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");
            CassonGavinCoeffs_.lookup("eta") >> eta_;                       // ADDED
            CassonGavinCoeffs_.lookup("tau_y") >> tau_y_;                   // ADDED
            CassonGavinCoeffs_.lookup("rho") >> rho_;                       // ADDED


Add coefficients in header file:

::

    $ nano CassonGavin.H
    
        dimensionedScalar eta_;      // ADDED
        dimensionedScalar tau_y_;    // ADDED
        dimensionedScalar rho_;      // ADDED

Compile in the main folder:

::

    $ wclean libso
    $ wmake libso


In constant/transportProperties

::

    transportModel  CassonGavin;

    CassonGavinCoeffs
    {
        eta        eta [ 1 -1 -1 0 0 0 0 ] 4.86;
        tau_y     tau_y [ 1 -1 -2 0 0 0 0 ] 14.38;
        rho        rho [1  -3 0 0 0 0 0]  1200; 
    }

In system/controlDict

::

    libs
    (
            "libCassonGavin.so"
    );

::

    $ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the velocity profile:

::

    $ sample
    
