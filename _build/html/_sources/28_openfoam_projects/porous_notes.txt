Porous Medium Notes
===================

.. contents::
   :local:

Re-run massflow = 40 laminar case
---------------------------------

::

	divSchemes
	{
		div((nuEff*dev2(T(grad(U))))) Gauss linear;
	}

	decomposeParDict

		numberOfSubdomains 8;

::

	$ plus

::

	$ decomposePar

::

	$ mpirun -np 8 buoyantBoussinesqSimpleFoam -parallel > log &

::

	$ tail -f log

Added switch in controlDict

::

	OptimisationSwitches
	{
		// Force dumping (at next timestep) and clean exit upon signal
		stopAtWriteNowSignal        10;
	}

::

	reconstructPar



File called plot.sh:

::

	#!/bin/bash

	foamLog log >/dev/null

	gnuplot -persist > /dev/null 2>&1 << EOF
			set logscale y
		    set format y "10^%+3T"
			set title "Residual vs. Iteration"
			set xlabel "Iteration"
			set ylabel "Residual"
			plot "logs/UxFinalRes_0" title 'Ux' with lines,\
			     "logs/UyFinalRes_0" title 'Uy' with lines,\
		         "logs/UzFinalRes_0" title 'Uz' with lines,\
		         "logs/TFinalRes_0" title 'T' with lines,\
		         "logs/p_rghFinalRes_1" title 'p rgh' with lines
	EOF


::

	chmod +x plot.sh

How to monitor residuals

::

	set logscale y
	set title "Residuals"
	set ylabel 'Residual'
	set xlabel 'Iteration'
	plot "< cat log | grep 'Solving for Ux' | cut -d' ' -f9 | tr -d ','" title 'Ux' with lines,\
	"< cat log | grep 'Solving for Uy' | cut -d' ' -f9 | tr -d ','" title 'Uy' with lines,\
	"< cat log | grep 'Solving for Uz' | cut -d' ' -f9 | tr -d ','" title 'Uz' with lines,\
	"< cat log | grep 'Solving for omega' | cut -d' ' -f9 | tr -d ','" title 'omega' with lines,\
	"< cat log | grep 'Solving for k' | cut -d' ' -f9 | tr -d ','" title 'k' with lines,\
	"< cat log | grep 'Solving for p' | cut -d' ' -f9 | tr -d ','" title 'p' with lines
	pause 1
	reread

::

	functions
	{
		  #include "pressureDifferencePatch"
	}

Command:

::

	gnuplot monitor -

How to kill mpirun

::

	kill -s 10 $MPIRUN_PID

Remove files:

::

	rm -rf 0.* [1-9]* postProcessing processor* *.foam logs

How to monitor pressure difference:

copy pressureDifferencePatch to system folder (same version as rest of code)
change <patch1> to inlet and <patch2> to outlet
chmod +x pressureDifferencePatch

Add to controlDict:

::

    functions
    {
        #include "pressureDifferencePatch"
    }

Command:

::

	gnuplot monitor2 -

STL File for CT Scan
--------------------

::

	$ grep 'solid'

Go to last line in nano

ctrl + w + v

Go to first line in nano

ctrl + w + y


MeshLab

Import > STL File

Unify duplicate vertices

File > Export mesh as > DXF


STL File
--------

Installed FreeCAD
Imported STL file
Part > Create Part from mesh



Tutorial on chtMultiRegionSimpleFoam
------------------------------------

Copied the version from OpenFOAM 2.3:

::

    http://openfoamwiki.net/index.php/Getting_started_with_chtMultiRegionSimpleFoam_-_planeWall2D

cd cht_001_2D_plate

$ ./Allrun

Boundaries may be of different types:
-patch (generic type)
-wall (for solid wall condition, useful for turbulence)
-cyclic (for cyclic simulations)
-symmetryPlane (for symmetry plane)
-empty (to specify that the simulation will be 2D or 1D)
-wedge (for axi-symmetric simulations)
-processor (for parallel computation, automatically defined during the decomposition domain process)


Generate the grid mesh: 
$ blockMesh
Check the mesh quality: 
$ checkMesh
View the mesh : 
$ paraFoam

Be careful with the units! In OpenFOAM incompressible solvers, the solved pressure is p = p'/rho

Plotting velocity vectors FROM CELL CENTRES (after running icoFoam) 
- conflict is that Paraview plots the vector at the face center, whereas OpenFOAM caluclates cell centre values

$ paraFoam -touch
$ paraFoam -builtin

Apply
Move to end of time
Properties > Representation > Surface
Properties > Coloring > U

Properties > Representation > Wireframe
Properties > Coloring > Solid Color

Switch to 2D

Filters > Search > Cell Centres > Apply

Filters > Search > Glyph > Apply
Properties > Glyph Type > 2D Glyph
Properties > Coloring > GlyphVector

Where are the solvers in OpenFOAM?

::

    $ cd $FOAM_APP/solvers
    $ ls

1. Darcy solver
---------------

Objective: develop a program that solves the velocity and pressure in a fully
saturated porous medium using Darcy's law.
div(U) = 0 (1)
U = -k/mu grad(p) (2)
How to solve this mathematical problem? The diffusion equation for the pressure field
is obtained by combining equation (1) and (2):

div(k/mu grad(p)) = 0

This equation is closed to the heat diffusion equation. Hence, we are going to program
our own solver on the basis of the existing laplacianFoam. To do so, we copy
laplacianFoam in our workspace.

::

    $ cd $WM_PROJECT_USER_DIR
    $ mkdir -p applications/solvers
    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $FOAM_APP/solvers/basic/laplacianFoam darcyFoam

Once the laplacianFoam solver has been copied into the user directory, we rename the main file and edit the Make/files    
    
::
    
    $ cd darcyFoam
    $ mv laplacianFoam.C darcyFoam.C
    $ gedit Make/files

Change files    
    
::

    darcyFoam.C

    EXE = $(FOAM_USER_APPBIN)/darcyFoam

We can now clean the previous compilation with wclean and compile this new program with wmake.

::

    $ wclean
    $ wmake
    
Edit the createFields.H file    
    
::

    $ gedit createFields.H

Edit createFields.H    
    
::

    Info<< "Reading field T\n" << endl;

    // Declaration of the pressure field p
    // * It is an instance of the object volScalarField (scalar field defined at the cells center),
    // * The file «p» must be read at the frist time step to satisfy the constructor. The initial values and boundary conditions are defined during the loading of 0/p.


    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE		
        ),
        mesh
    );

    // Declaration of the velocity vector field U
    // * It is an instance of the object volVectorField (vector field defined at the cells center),
    // * U is not read from a file (even if 0/U exist)
    // * To satisfy the constructor of the object volVectorField, units and initial values are defined with an additional argument. By default, the boundary conditions are zeroGradient,
    // * The file « U » will be written at every output times.

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),vector::zero)
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Declaration of the fluid viscosity mu and the permeability k.
    // They will be loaded from « constant/transportProperties »

    Info<< "Reading diffusivity k\n" << endl;

    dimensionedScalar k
    (
        transportProperties.lookup("k")
    );

    Info<< "Reading fluid viscosity mu\n" << endl;

    dimensionedScalar mu
    (
        transportProperties.lookup("mu")
    );
    
::

    $ gedit darcyFoam.C    
    
    
::

    #include "fvCFD.H"
    #include "simpleControl.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    int main(int argc, char *argv[])
    {
        #include "setRootCase.H"

        #include "createTime.H"
        #include "createMesh.H"
        #include "createFields.H"

        simpleControl simple(mesh);

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        Info<< "\nCalculating temperature distribution\n" << endl;

        while (simple.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            while (simple.correctNonOrthogonal())
            {

                // The pressure field p is solved implicitly by a diffusion equation
                solve
                (
                    fvm::laplacian(k/mu, p)
                );
            }
            
            // The velocity vector U is deduced from the pressure field using the Darcy's law.
            U = -k/mu*fvc::grad(p);

            // runtime selectable variables
            runTime.write();

            // This file contained references to temperature gradient, which is not needed
            //#include "write.H" 

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }

        Info<< "End\n" << endl;

        return 0;
    }
    

The useless files are removed and the darcyFoam executable is then compiled.

::

    $ rm write.H
    $ wclean
    $ wmake

To prepare this « case » we are going to use the tutorial laplacianFoam/flange. The setup will be relatively similar since this latter also solves a diffusion equation.

::

    $ run
    $ cp -r $FOAM_TUTORIALS/basic/laplacianFoam/flange .
    $ mv flange porous_005_darcyFoam
    $ cd porous_005_darcyFoam
    $ rm Allrun Allclean flange.ans
    
To save time, we can pick up and modify an existing blockMeshDict

::

$ cp -rf $FOAM_TUTORIALS/incompressible/icoFoam/cavity/constant/polyMesh/blockMeshDict constant/polyMesh/blockMeshDict
$ gedit constant/polyMesh/blockMeshDict
    

::    
    
    convertToMeters 1;

    // The vertices can be defined using variables Lx, Ly, Lz. It saves time to modify the size of the domain

    Lx 10;
    Ly 0.1;
    Lz 0.1;

    vertices
    (
        (0 0 0)
        ($Lx 0 0)
        ($Lx $Ly 0)
        (0 $Ly 0)
        (0 0 $Lz)
        ($Lx 0 $Lz)
        ($Lx $Ly $Lz)
        (0 $Ly $Lz)
    );

    // Mesh definition (homogeneous grid with a single cell in the y and z axis since the simulation is 1D)

    blocks
    (
        hex (0 1 2 3 4 5 6 7) (60 1 1) simpleGrading (1 1 1)
    );

    edges
    (
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

        frontAndBack
        {
            type empty;
            faces
            (
                            (1 5 4 0)
                            (3 7 6 2)
                (0 3 2 1)
                (4 5 6 7)
            );
        }
    );

    mergePatchPairs
    (
    );
   
The grid is generated using blockMesh

::

    $ blockMesh    

The initial and boundary conditions are changed:

::

    $ mv 0/T 0/p
    $ gedit 0/p   

Change the temperature file to a pressure file    
    
::

    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volScalarField;
        object      p;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [1 -1 -2 0 0 0 0];

    internalField   uniform 0;

    boundaryField
    {
    
        // A pressure drop is imposed between the inlet and the outlet of the computational domain
        inlet
        {
            type            fixedValue;
            value	    uniform 1e2;
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
    
Change the transportProperties:

::

    $ gedit constant/transportProperties

::

    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "constant";
        object      transportProperties;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    k              k  [ 0  2  0 0 0 0 0 ] 1e-09;
    mu             mu [ 1 -1 -1 0 0 0 0 ] 1e-05;
    
Change the control dictionary:

::

    $ gedit system/controlDict

::

    // Since darcyFoam is a steady-state solver without relaxation factor, only one time step is necessary.

    endTime         1;

    deltaT          1;

Change the fvSchemes for the laplacian of k, mu and p.

Remove the fluxRequired for temperature.

::

    laplacianSchemes
    {
        default         none;
        laplacian((k|mu),p) Gauss linear corrected;
    }


    //fluxRequired
    //{
    //    default         no;
    //    T               ;
    //}

Solve for pressure not temperature    

::

    solvers
    {
        p
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-06;
            relTol          0;
        }
    }


Run the simulation:

::

    $ darcyFoam

Results will be plotted using the sample utility and the program Gnuplot. As blockMesh, the program sample requires an input dictionary located in /system    
    
::
    
    $ cp $FOAM_UTILITIES/postProcessing/sampling/sample/sampleDict system/.
    $ gedit system/sampleDict    

Edit the sampleDict    
    
::

    sets
    (
        lineX1
        {
            type        midPoint;
            axis        distance;

            //- cavity. Slightly perturbed so not to align with face or edge.
            start       (0 0.05 0.005);
            end         (10 0.05 0.005);
        }

    );
    
Run the sample dictionary (it won't be able to find 0/U):

::

    $ sample
    

Plot the results    
    
::

    $ gnuplot
    > set xlabel "distance (m)"
    > set ylabel "Pressure (kg/m/s)"
    > plot "postProcessing/sets/1/lineX1_p.xy" using 1:2 with lines lw 4 title "p"
    
2. Darcy solver with heat transfer
----------------------------------

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r darcyFoam darcyFoamTemperature
    $ cd darcyFoamTemperature
    $ mv darcyFoam.C darcyFoamTemperature.C
    $ gedit Make/files
    
Changes files file    
    
::    
    
    darcyFoamTemperature.C

    EXE = $(FOAM_USER_APPBIN)/darcyFoamTemperature    

::

    $ wclean
    $ wmake
    
::

    $ gedit createFields.H

Added phi and temperature    
    
::

    // Declaration of the velocity flux phi.
    // * It is a surface field (U is projected onto the face of each cell of the grid)
    // * It is necessary when using the divergence opereator (fvm::div(phi,T) )
    // * Can also be declared using #include "createPhi.H"


        surfaceScalarField phi("phi", linearInterpolate(U) & mesh.Sf());

        Info<< "Reading field T\n" << endl;

    // Declaration of the temperature field T . (Do not copy everything manually: copy/paste the declaration of volScalarField p and replace p by T)


        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE		
            ),
            mesh
        );

Declare thermal conductivity, porosity, heat capacities of fluid and solid.

::

    // Beside the viscosity mu and the permeability k of the porous medium, we also declare the thermal conductivity DT, the porosity eps and the heat capacities rhoCps and rhoCpf. They are loaded from the file « constant/transportProperties »

    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

    Info<< "Reading fluid viscosity eps\n" << endl;

    dimensionedScalar eps
    (
        transportProperties.lookup("eps")
    );

    Info<< "Reading fluid viscosity rhoCps\n" << endl;

    dimensionedScalar rhoCps
    (
        transportProperties.lookup("rhoCps")
    );
    
    Info<< "Reading fluid viscosity rhoCpf\n" << endl;

    dimensionedScalar rhoCpf
    (
        transportProperties.lookup("rhoCpf")
    );   
    
Re-make the solver

::

    $ wclean
    $ wmake
    
Add the temperature equation:

::

    // The surface flux phi is updated from the new value of the velocity profile U.
    phi = linearInterpolate(U) & mesh.Sf();

    // Solve the advection/diffusion equation for the temperature transport
    solve
    (
            (eps*rhoCpf+(1.0-eps)*rhoCps)*fvm::ddt(T) 
            + rhoCpf*fvm::div(phi,T)
            ==	
            fvm::laplacian(DT,T)
    );

Re-make the solver

::

    $ wclean
    $ wmake
    
    
Copy the previous exercise

::

    $ run
    $ cp -r porous_005_darcyFoam porous_006_darcyFoamTemperature
    $ cd porous_006_darcyFoamTemperature
    $ rm -rf 1 postProcessing *.foam
    
Create temperature boundary conditions:

::

    $ cp 0/p 0/T
    $ gedit 0/T
    
::

    dimensions      [0 0 0 1 0 0 0];

    internalField   uniform 273;

    boundaryField
    {
        inlet
        {
            type            fixedValue;
            value	    uniform 573;
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
    
Edit the transport properties

::

    $ gedit constant/transportProperties
    
    
::

    k              k  [ 0  2  0 0 0 0 0 ] 1e-09;

    mu             mu [ 1 -1 -1 0 0 0 0 ] 1e-05;

    eps            eps [ 0 0 0 0 0 0 0 ] 0.4;

    DT             DT [ 1 1 -3 -1 0 0 0 ] 1e-02;

    rhoCps         rhoCps [ 1 -1 -2 -1 0 0 0 ] 2e4;

    rhoCpf         rhoCpf [ 1 -1 -2 -1 0 0 0 ] 5e3;
    
    
Edit the solver properties add temperature:

::

    solvers
    {
        p
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-06;
            relTol          0;
        }

        T
        {
            solver          PBiCG;
            preconditioner  DILU;
            tolerance       1e-06;
            relTol          0;
        }
    }

    SIMPLE
    {
        nNonOrthogonalCorrectors 2;
    }
    
Edit the discretisation fvSchemes
    
::

    ddtSchemes
    {
        default         Euler;
    }

    gradSchemes
    {
        default         Gauss linear;
        grad(T)         Gauss linear;
    }

    divSchemes
    {
        default         none;
            div(phi,T)		Gauss linear;
    }

    laplacianSchemes
    {
        default         none;
        laplacian((k|mu),p) Gauss linear corrected;
            laplacian(DT,T) Gauss linear corrected;
    }

    interpolationSchemes
    {
        default         linear;
    }

    snGradSchemes
    {
        default         corrected;
    }

    fluxRequired
    {
        default         no;
        T               ;
    }    
    
Set a probe to record temperature versus time

::

    $ gedit system/controlDict
    

The probes are functions that are executed on-the-fly during the simulation.
They allow to record the temperature evolution vs time at the probe location.
You can specify as many probes as you want.        
    
::

    // darcyFoamTemperature is a temporal solver

    endTime         60000;

    deltaT          100;

    writeControl    runTime;

    writeInterval   1000;


    functions
    {
        probes
        {
            type 		probes;
            functionObjectLibs 	("libsampling.so");
            enabled 		true;
            outputControl 	timeStep;
            outputInterval 	1;

            fields
            (
                    T
            );
            
            probeLocations
            (
                    (2 0.05 0.05) // Probe 1
                    (5 0.05 0.05) // Probe 2
                    (9 0.05 0.05) // Probe 3
            );
        }
    }
    
Plot the temperature, create a new file, plot_probes

::

    set key at 50000, 400
    set ylabel "temperature (K)"
    set xlabel "time (s)"

    plot "postProcessing/probes/0/T" using 1:2 with lines lw 4 title "Probe x=2m" ,\
    "postProcessing/probes/0/T" using 1:3 with lines lw 4 title "Probe x=5m" ,\
    "postProcessing/probes/0/T" using 1:4 with lines lw 4 title "Probe x=9m"

Run the script:

::

    $ gnuplot -persist plot_probes
 
Changed scheme to vanLeer for convective term - gives smoother curves

3. Customize boundary conditions
--------------------------------

Where is the source code for boundary conditions?

::

    $ cd $FOAM_SRC/finiteVolume/fields/fvPatchFields
    $ ls
    
How are boundary conditions defined?

To define boundary conditions that depends on time or on other variables, there are several
possibilities:
- Hardcoding in the solver,
- Program customized boundary conditions,
- With an additional package such as swak4Foam
    
The porous medium was solved for the pressure, so the boundary conditions for the preussure must be specified. 

However, it is more convienient to specify the velocity, hence n.grad(p) = -(mu/k) n.U

Create a new condition inspired by fixedFluxPressure

::

    $ mkdir -p $WM_PROJECT_USER_DIR/boundary/
    $ cd $WM_PROJECT_USER_DIR/boundary/
    $ cp -r $FOAM_SRC/finiteVolume/fields/fvPatchFields/derived/fixedFluxPressure darcyGradPressure
    $ cd darcyGradPressure

Replace strings fixedFluxPressure with darcyGradPressure
    
::

    $ rename 's/fixedFluxPressure/darcyGradPressure/g' *.*
    $ sed -i 's/fixedFluxPressure/darcyGradPressure/g' *.*
    
Create a files file in Make directory

::

    $ mkdir Make
    $ gedit Make/files
    
::
    
    darcyGradPressureFvPatchScalarField.C

    LIB = $(FOAM_USER_LIBBIN)/ldarcyGradPressure
       
    
::

    $ gedit Make/options    
    
::
    
    EXE_INC = -I$(LIB_SRC)/finiteVolume/lnInclude

    LIB_LIBS = -lfiniteVolume

Clean and make the files

::

    $ wclean
    $ wmake

Declare `UName_` in darcyGradPressureFvPatchScalarField.H

::

	/*---------------------------------------------------------------------------*\
		         Class darcyGradPressureFvPatchScalarField Declaration
	\*---------------------------------------------------------------------------*/

	class darcyGradPressureFvPatchScalarField
	:
		public fixedGradientFvPatchScalarField      // This boundary condition derives from the class fixedGradient
	{
		// Private data

		    //- Name of the velocity field to calculate grad(p)
		    word UName_;

Add UName in darcyGradPressureFvPatchScalarField.C:

::

	#include "darcyGradPressureFvPatchScalarField.H"
	#include "fvPatchFieldMapper.H"
	#include "volFields.H"
	#include "surfaceFields.H"
	#include "addToRunTimeSelectionTable.H"

	// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

	Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		UName_("U")		// ADDED
	{}


	Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
	(
		const darcyGradPressureFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		UName_(ptf.UName_)		// ADDED
	{
		// REMOVED:
		//patchType() = ptf.patchType();

		// Map gradient. Set unmapped values and overwrite with mapped ptf
		//gradient() = 0.0;
		//gradient().map(ptf.gradient(), mapper);

		// Evaluate the value field from the gradient if the internal field is valid
		/*
		if (&iF && iF.size())
		{
		    scalarField::operator=
		    (
		        //patchInternalField() + gradient()/patch().deltaCoeffs()
		        // ***HGW Hack to avoid the construction of mesh.deltaCoeffs
		        // which fails for AMI patches for some mapping operations
		        patchInternalField() + gradient()*(patch().nf() & patch().delta())
		    );
		}
		*/
	}


	Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		UName_(dict.lookupOrDefault<word>("U","U"))			// ADDED
	{
		// REMOVED:

		/*
		if (dict.found("value") && dict.found("gradient"))
		{
		    fvPatchField<scalar>::operator=
		    (
		        scalarField("value", dict, p.size())
		    );
		    gradient() = scalarField("gradient", dict, p.size());
		}
		else
		{
		*/
		    fvPatchField<scalar>::operator=(patchInternalField());
		    gradient() = 0.0;
		//}
	}


	Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
	(
		const darcyGradPressureFvPatchScalarField& wbppsf
	)
	:
		fixedGradientFvPatchScalarField(wbppsf),
		UName_(wbppsf.UName_)		// ADDED
	{}


	Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
	(
		const darcyGradPressureFvPatchScalarField& wbppsf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(wbppsf, iF),
		UName_(wbppsf.UName_)		// ADDED
	{}




Change the evaluation of the gradient:

::

	// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

	void Foam::darcyGradPressureFvPatchScalarField::updateCoeffs
	(
		const scalarField& snGradp
	)
	{
		if (updated())
		{
		    return;
		}

		// We recover the value of U at the boundary:

		const fvPatchField<vector>& U = 
			patch().lookupPatchField<volVectorField, vector>(UName_);

		// Extract the dictionary from the database:
		const dictionary& transportProperties = db().lookupObject<IOdictionary> ("transportProperties");

		// Extract viscosity and permeability:
		dimensionedScalar mu(transportProperties.lookup("mu"));
		dimensionedScalar k(transportProperties.lookup("k"));
	
		// The pressure gradient is evaluated at the boundary with the formula:
		// n.grad(p) = -(mu/k)n.U
		// mu.value() allows the access of the value of the object mu declared as a dimensionedScalar
		// patch().nf() returns the normal vector to the patch
		gradient() = -mu.value()/k.value()*(U & patch().nf());

	 	// curTimeIndex_ = this->db().time().timeIndex();

		// gradient() = snGradp;
		fixedGradientFvPatchScalarField::updateCoeffs();
	}


	void Foam::darcyGradPressureFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
		    return;
		}

		// REMOVED:
		/*
		if (curTimeIndex_ != this->db().time().timeIndex())
		{
		    FatalErrorIn("darcyGradPressureFvPatchScalarField::updateCoeffs()")
		        << "updateCoeffs(const scalarField& snGradp) MUST be called before"
		           " updateCoeffs() or evaluate() to set the boundary gradient."
		        << exit(FatalError);
		}
		*/
	}

The library ldarcyGradPressure.so is now compiled and available for all the solvers

::

    $ wclean
    $ wmake


4. Two-equation model
---------------------

This solver will be based on darcyFoamTemperature

::

	$ cd $WM_PROJECT_USER_DIR/applications/solvers/
	$ mkdir darcyFoamTemperatureTwo
	$ cp -r darcyFoamTemperature darcyFoamTemperatureTwo
	$ cd darcyFoamTemperatureTwo
	$ mv darcyFoamTemperature.C darcyFoamTemperatureTwo.C

Edit the path for the executable:

::

	$ gedit Make/files

Change the files file:

::

	darcyFoamTemperatureTwo.C

	EXE = $(FOAM_USER_APPBIN)/darcyFoamTemperatureTwo

Create the fields U, Ts and Tf:

::

	$ gedit createFields.H

::

    Info<< "Reading field p\n" << endl;

    // Declaration of the pressure field p
    // * It is an instance of the object volScalarField (scalar field defined at the cells center),
    // * The file «p» must be read at the frist time step to satisfy the constructor. The initial values and boundary conditions are defined during the loading of 0/p.


            volScalarField p
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE		
                ),
                mesh
            );

    // Declaration of the velocity vector field U
    // * It is an instance of the object volVectorField (vector field defined at the cells center),
    // * U is not read from a file (even if 0/U exist)
    // * To satisfy the constructor of the object volVectorField, units and initial values are defined with an additional argument. By default, the boundary conditions are zeroGradient,
    // * The file « U » will be written at every output times.

            Info<< "Reading field U\n" << endl;

    // The field U is now initialised from 0/U, which allows us to define the boundary condition for U

            volVectorField U
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,			 
                    IOobject::AUTO_WRITE
                ),
                mesh,
                    ///dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0),vector::zero)
            );

    // phi is created by calling the createPhi.H

    #include "createPhi.H"



    // Declaration of the velocity flux phi.
    // * It is a surface field (U is projected onto the face of each cell of the grid)
    // * It is necessary when using the divergence opereator (fvm::div(phi,T) )
    // * Can also be declared using #include "createPhi.H"


    //	surfaceScalarField phi("phi", linearInterpolate(U) & mesh.Sf()); // REMOVED

            // Declaration of the temperature field Ts 

            Info<< "Reading field Ts\n" << endl;

            volScalarField Ts
            (
                IOobject
                (
                    "Ts",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE		
                ),
                mesh
            );

    // Declaration of the temperature field Tf 

            Info<< "Reading field Tf\n" << endl;

            volScalarField Tf
            (
                IOobject
                (
                    "Tf",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE		
                ),
                mesh
            );

            Info<< "Reading transportProperties\n" << endl;

            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            );

    // Beside the viscosity mu and the permeability k of the porous medium, we also declare the thermal conductivity DT, the porosity eps and the heat capacities rhoCps and rhoCpf. They are loaded from the file « constant/transportProperties »

            Info<< "Reading permeability k\n" << endl;

            dimensionedScalar k
            (
                transportProperties.lookup("k")
            );

            Info<< "Reading permeability mu\n" << endl;

            dimensionedScalar mu
            (
                transportProperties.lookup("mu")
            );

            Info<< "Reading diffusivity DTf\n" << endl;

            dimensionedScalar DTf
            (
                transportProperties.lookup("DTf")
            );

            Info<< "Reading diffusivity DTs\n" << endl;

            dimensionedScalar DTs
            (
                transportProperties.lookup("DTs")
            );

            Info<< "Reading fluid viscosity eps\n" << endl;

            dimensionedScalar eps
            (
                transportProperties.lookup("eps")
            );

            Info<< "Reading fluid viscosity rhoCps\n" << endl;

            dimensionedScalar rhoCps
            (
                transportProperties.lookup("rhoCps")
            );
            
            Info<< "Reading fluid viscosity rhoCpf\n" << endl;

            dimensionedScalar rhoCpf
            (
                transportProperties.lookup("rhoCpf")
            );  

            Info<< "Reading heat exchange coefficient h\n" << endl;

            dimensionedScalar h
            (
                transportProperties.lookup("h")
            ); 

    // Declaration of the fluid viscosity mu and the permeability k.
    // They will be loaded from « constant/transportProperties »
    /*
            Info<< "Reading diffusivity k\n" << endl;

            dimensionedScalar k
            (
                transportProperties.lookup("k")
            );

            Info<< "Reading fluid viscosity mu\n" << endl;

            dimensionedScalar mu
            (
                transportProperties.lookup("mu")
            );
    */

Add the two-equations to darcyFoamTemperatureTwo.C

Correct boundary conditions

::

	$ gedit darcyFoamTemperatureTwo.C

::

	#include "fvCFD.H"
	#include "simpleControl.H"

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	int main(int argc, char *argv[])
	{
		#include "setRootCase.H"

		#include "createTime.H"
		#include "createMesh.H"
		#include "createFields.H"

		simpleControl simple(mesh);

		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

		Info<< "\nCalculating temperature distribution\n" << endl;

		while (simple.loop())
		{
		    Info<< "Time = " << runTime.timeName() << nl << endl;

		    while (simple.correctNonOrthogonal())
		    {

				// The pressure field p is solved implicitly by a diffusion equation
				// The matrix for the pressure pEqn is decalred as a fvScalarMatrix and contructed from the discretisation with the finite volume method of the laplacian operator. 
				// This type of dicretisation is more flexible than solve(fvm::laplacian(...))
				fvScalarMatrix pEqn
		        (
		            fvm::laplacian(k/mu, p)
		        );

				// The matrix is inversed with pEqn.solve()
				pEqn.solve();

				if(simple.finalNonOrthogonalIter())
				{
					phi = -pEqn.flux();
				}
		    }


			// The velocity vector U is deduced from the pressure field using the Darcy's law.
			U = -k/mu*fvc::grad(p);
			// The boundary conditions may have been altered and do not correspond to what is specified in 0/U
			// This function means the BCs are those specified in 0/U
			U.correctBoundaryConditions();

			// The surface flux phi is updated from the new value of the velocity profile U.
			//phi = linearInterpolate(U) & mesh.Sf();

			// Solve the temperature in the fluid

			// The laplacian discretization scheme is specified in system/fvSchemes in front of the keyword « laplacian(DT,T) ». A part of the exchange term is treated implicitly, the other part explicitly.
			fvScalarMatrix TfEqn
			(
				eps*rhoCpf*fvm::ddt(Tf) + rhoCpf*fvm::div(phi,Tf)
				==
				fvm::laplacian(eps*DTf,Tf,"laplacian(DT,T)") - fvm::Sp(h,Tf) + h*Ts
			);
			TfEqn.solve();

			// Solve the temperature in the solid
			fvScalarMatrix TsEqn
			(
				(1.0-eps)*rhoCps*fvm::ddt(Ts)
				==
				fvm::laplacian((1.0-eps)*DTs,Ts,"laplacian(DT,T)") - fvm::Sp(h,Ts) + h*Tf
			);
			TsEqn.solve();



			// runtime selectable variables
			runTime.write();

			// This file contained references to temperature gradient, which is not needed
		    //#include "write.H" 

		    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
		        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
		        << nl << endl;
		}

		Info<< "End\n" << endl;

		return 0;
	}

5. Run the Two Temperature Case
-------------------------------

Copy the previous case:

::

    $ run
    $ cp -r porous_007_darcyFoamTemperatureVanLeer porous_008_darcyFoamTemperatureTwo
    $ cd porous_008_darcyFoamTemperatureTwo

Remove folders:

::

    $ rm -rf 0.* [1-9]* postProcessing

Create new initial conditions:

::

	$ mv 0/T 0/Tf
	$ cp 0/Tf 0/Ts
	$ cp 0/p 0/U

Edit the 0/U file:

::

	FoamFile
	{
		version     2.0;
		format      ascii;
		class       volVectorField;
		object      U;
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	dimensions      [0 1 -1 0 0 0 0];

	internalField   uniform (0 0 0);

	boundaryField
	{
		inlet
		{
		    type            fixedValue;
			value			uniform (1.4e-4 0 0);
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

	// ************************************************************************* //

Edit the 0/p file - test with 10Pa inlet pressure:

::

	FoamFile
	{
		version     2.0;
		format      ascii;
		class       volScalarField;
		object      p;
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	dimensions      [1 -1 -2 0 0 0 0];

	internalField   uniform 0;

	boundaryField
	{
		inlet
		{
		    type            fixedValue;
			//type			darcyGradPressure;
			value			uniform 10;
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

	// ************************************************************************* //

Edit the 0/Ts file:

::

	FoamFile
	{
		version     2.0;
		format      ascii;
		class       volScalarField;
		object      Ts;
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	dimensions      [0 0 0 1 0 0 0];

	internalField   uniform 573;

	boundaryField
	{
		inlet
		{
		    type        zeroGradient;
		    //value	    uniform 573;
		}

		outlet
		{
		    type        zeroGradient;
		}

		frontAndBack
		{
		    type        empty;
		}
	}

	// ************************************************************************* //

Edit the 0/Tf file:

::

	FoamFile
	{
		version     2.0;
		format      ascii;
		class       volScalarField;
		object      Tf;
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	dimensions      [0 0 0 1 0 0 0];

	internalField   uniform 573;

	boundaryField
	{
		inlet
		{
		    type        fixedValue;
		    value	    uniform 273;
		}

		outlet
		{
		    type        zeroGradient;
		}

		frontAndBack
		{
		    type        empty;
		}
	}

	// ************************************************************************* //


Edit the transport properties

::

	$ gedit constant/transportProperties


::

	FoamFile
	{
		version     2.0;
		format      ascii;
		class       dictionary;
		location    "constant";
		object      transportProperties;
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	k              	k  		[ 0  2  0  0 0 0 0 ] 1e-09;

	mu             	mu 		[ 1 -1 -1  0 0 0 0 ] 1e-05;

	eps            	eps 	[ 0  0  0  0 0 0 0 ] 0.4;

	DTf             DTf 	[ 1  1 -3 -1 0 0 0 ] 1e-04;

	DTs             DTs 	[ 1  1 -3 -1 0 0 0 ] 1e-02;

	rhoCps         	rhoCps 	[ 1 -1 -2 -1 0 0 0 ] 2e4;

	rhoCpf         	rhoCpf 	[ 1 -1 -2 -1 0 0 0 ] 5e3;

	h				h		[ 1 -1 -3 -1 0 0 0 ] 5e-1;


	// ************************************************************************* //


Change the solvers in fvSolution:

::

	solvers
	{
		p
		{
		    solver          PCG;
		    preconditioner  DIC;
		    tolerance       1e-06;
		    relTol          0;
		}

		Ts
		{
		    solver          PCG;
		    preconditioner  DIC;
		    tolerance       1e-06;
		    relTol          0;
		}

		Tf
		{
		    solver          PBiCG;
		    preconditioner  DILU;
		    tolerance       1e-06;
		    relTol          0;
		}
	}

	SIMPLE
	{
		nNonOrthogonalCorrectors 2;
	}    
		
Change the schemes in fvSchemes:

::

	FoamFile
	{
		version     2.0;
		format      ascii;
		class       dictionary;
		location    "system";
		object      fvSchemes;
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	ddtSchemes
	{
		default         Euler;
	}

	gradSchemes
	{
		default         Gauss linear;
		grad(T)         Gauss linear;
	}

	divSchemes
	{
		default         none;
		div(phi,Tf)		Gauss vanLeer;
	}

	laplacianSchemes
	{
		default         	none;
		laplacian((k|mu),p) Gauss linear corrected;
		laplacian(DT,T) 	Gauss linear corrected;
	}

	interpolationSchemes
	{
		default         linear;
	}

	snGradSchemes
	{
		default         corrected;
	}

	fluxRequired
	{
		p;
	}


Change sampleDict to sample fluid and solid temperatures
    
::

	fields
	(
		Tf
		Ts
	);

	sets
	(
		lineX1
		{
		    type        midPoint;
		    axis        distance;

		    //- cavity. Slightly perturbed so not to align with face or edge.
		    start       (0 0.05 0.005);
		    end         (10 0.05 0.005);
		}

	);


Delete the old files and re-run:

::

	$ rm -rf 0.* [1-9]* postProcessing
	$ darcyFoamTemperatureTwo
	$ sample
	$ gnuplot -persist plot_probes


5. Test with New Boundary Condition
-----------------------------------

Copy the previous case:

::

    $ run
    $ cp -r porous_008_darcyFoamTemperatureTwo porous_009_darcyFoamTemperatureTwoBC
    $ cd porous_009_darcyFoamTemperatureTwoBC

Edit the pressure file to include the new condition:

::

	$ gedit 0/p

::

	FoamFile
	{
		version     2.0;
		format      ascii;
		class       volScalarField;
		object      p;
	}
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	dimensions      [1 -1 -2 0 0 0 0];

	internalField   uniform 0;

	boundaryField
	{
		inlet
		{
		    //type            fixedValue;
			type			darcyGradPressure;
			value			uniform 0;
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

	// ************************************************************************* //

Add the library to controlDict:

::

    libs ("ldarcyGradPressure.so");

    functions
    {
        probes
        {
            type 		probes;
            functionObjectLibs 	("libsampling.so");
            enabled 		true;
            outputControl 	timeStep;
            outputInterval 	1;

            fields
            (
                    T
            );
            
            probeLocations
            (
                    (2 0.05 0.05) // Probe 1
                    (5 0.05 0.05) // Probe 2
                    (9 0.05 0.05) // Probe 3
            );
        }
    }

Delete the old files and re-run:

::

	$ rm -rf 0.* [1-9]* postProcessing
	$ darcyFoamTemperatureTwo
	$ sample
	$ gnuplot -persist plot_probes


XXX. Make porousSimpleFoamCopy
------------------------------

Simple Foam

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $FOAM_APP/solvers/incompressible/simpleFoam simpleCopyFoam
    
    $ cd simpleCopyFoam
    $ nano Make/files
        simpleCopyFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/simpleCopyFoam
    $ mv simpleFoam.C simpleCopyFoam.C
    $ wclean
    $ wmake

Porous Simple Foam

::

    $ rm -rf SRFSimpleFoam
    $ mv porousSimpleFoam porousSimpleCopyFoam
    $ cd porousSimpleCopyFoam
    $ nano Make/files
    
        porousSimpleCopyFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/porousSimpleCopyFoam
    
    $ mv porousSimpleFoam.C porousSimpleCopyFoam.C
    $ wclean
    $ wmake
    
Setup a 1D Case

::

    $ run
    $ cp -rf porous_005_darcyFoam porous_010_porousSimpleCopyFoam
    $ rm -rf 0.* [1-9]* postProcessing

Add a zone for porosity
    
::

    $ nano constant/polyMesh/blockMeshDict
    
        blocks
        (
            hex (0 1 2 3 4 5 6 7) 
            porosity (60 1 1) simpleGrading (1 1 1)
        );

Change solver
    
::
    
    $ nano system/controlDict
        application porousSimpleCopyFoam
    $ nano 0/p
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

    $ cp 0/p 0/U
    $ nano 0/U

        FoamFile
        {
            version     2.0;
            format      ascii;
            class       volVectorField;
            object      U;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        dimensions      [0 1 -1 0 0 0 0];

        internalField   uniform (0 0 0);

        boundaryField
        {
            inlet
            {
                type            fixedValue;
                value           uniform (1.4e-4 0 0);
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

    $ nano system/controlDict
        libs ("ldarcyGradPressure.so");

Go through each dictionary and match it with angleDuctExplicit (as this used porousSimpleFoam)

::

    $ nano constant/transportProperties
    
        transportModel  Newtonian;
        nu              nu [0 2 -1 0 0 0 0] 1.5e-05;

    $ nano constant/RASProperties
    
        RASModel        laminar;

        turbulence      off;

        printCoeffs     off;
    
    $ nano constant/porosityProperties
    
    
        porosity1
        {
            type            DarcyForchheimer;
            active          yes;
            cellZone        porosity;

            DarcyForchheimerCoeffs
            {
                d   d [0 -2 0 0 0 0 0] (1e9 0 0);
                f   f [0 -1 0 0 0 0 0] (0 0 0);

                coordinateSystem
                {
                    type    cartesian;
                    origin  (0 0 0);
                    coordinateRotation
                    {
                        type    axesRotation;
                        e1      (1 0 0);
                        e2      (0 1 0);
                    }
                }
            }
        }
    
    $ nano system/controlDict
        
        application     porousSimpleCopyFoam;

        startFrom       startTime;

        startTime       0;

        stopAt          endTime;

        endTime         200;

        deltaT          1;

        writeControl    timeStep;

        writeInterval   10;

        purgeWrite      0;

        writeFormat     ascii;

        writePrecision  6;

        writeCompression off;

        timeFormat      general;

        timePrecision   6;

        runTimeModifiable true;

        libs ("ldarcyGradPressure.so");
    
    $ nano system/fvSchemes    
    
        ddtSchemes
        {
            default         steadyState;
        }

        gradSchemes
        {
            default         Gauss linear;
        }

        divSchemes
        {
            default         none;

            div(phi,U)      bounded Gauss upwind;
            div((nuEff*dev(T(grad(U))))) Gauss linear;
            div(phi,epsilon) bounded Gauss upwind;
            div(phi,k)      bounded Gauss upwind;
        }

        laplacianSchemes
        {
            default         Gauss linear corrected;
        }

        interpolationSchemes
        {
            default         linear;
        }

        snGradSchemes
        {
            default         corrected;
        }

        fluxRequired
        {
            default         no;
            p               ;
        }
    
    $ nano system/fvSchemes     
    
        application     porousSimpleCopyFoam;

        startFrom       startTime;

        startTime       0;

        stopAt          endTime;

        endTime         200;

        deltaT          1;

        writeControl    timeStep;

        writeInterval   10;

        purgeWrite      0;

        writeFormat     ascii;

        writePrecision  6;

        writeCompression off;

        timeFormat      general;

        timePrecision   6;

        runTimeModifiable true;

        libs ("ldarcyGradPressure.so"); 

Test Darcy BC:
    
::

    $ cp -rf porous_010_porousSimpleCopyFoam porous_011_porousSimpleCopyDarcyBCFoam
    $ cd porous_011_porousSimpleCopyDarcyBCFoam    
    $ rm -rf 0.* [1-9]* postProcessing *.foam   
    
    $ nano 0/p
        
        inlet
        {
            type            darcyGradPressure;
            value           uniform 0;
        }
        
    $ constant/transportProperties
    
        // For Darcy BC:
        k               k  [0 2   0 0 0 0 0] 1e-09;
        mu              mu [1 -1 -1 0 0 0 0] 1e-05;

    $ porousSimpleCopyFoam
    $ sample
    $ gnuplot -persist plot_figure
    
    
Add temperature to solver (local equilibirum assumption):

Simple Foam

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $FOAM_APP/solvers/incompressible/simpleFoam simpleTemperatureFoam
    
    $ cd simpleTemperatureFoam
    $ nano Make/files
        simpleTemperatureFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/simpleTemperatureFoam
    $ mv simpleFoam.C simpleTemperatureFoam.C
    $ wclean
    $ wmake

Porous Simple Foam

::

    $ rm -rf SRFSimpleFoam
    $ mv porousSimpleFoam porousSimpleTemperatureFoam
    $ cd porousSimpleTemperatureFoam
    $ nano Make/files
    
        porousSimpleTemperatureFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/porousSimpleTemperatureFoam
    
    $ mv porousSimpleFoam.C porousSimpleTemperatureFoam.C
    $ wclean
    $ wmake

Add temperature field to createFields.H

::

	$ nano createFields.H

		Info<< "Reading field T\n" << endl;
		volScalarField T
		(
		    IOobject
		    (
		        "T",
		        runTime.timeName(),
		        mesh,
		        IOobject::MUST_READ,
		        IOobject::AUTO_WRITE
		    ),
		    mesh
		);

		#include "createPhi.H"

Add coefficients

::

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

    Info<< "Reading fluid viscosity eps\n" << endl;

    dimensionedScalar eps
    (
        transportProperties.lookup("eps")
    );

    Info<< "Reading fluid viscosity rhoCps\n" << endl;

    dimensionedScalar rhoCps
    (
        transportProperties.lookup("rhoCps")
    );
    
    Info<< "Reading fluid viscosity rhoCpf\n" << endl;

    dimensionedScalar rhoCpf
    (
        transportProperties.lookup("rhoCpf")
    );   

Add TEqn.H to porousSimpleTemperatureFoam.C

::

    turbulence->correct();

    #include "TEqn.H"

Create TEqn.H

::

	cp UEqn.H TEqn.H

Add this to TEqn.H

::

    solve
    (
            (eps*rhoCpf+(1.0-eps)*rhoCps)*fvm::ddt(T) 
            + rhoCpf*fvm::div(phi,T)
            ==	
            fvm::laplacian(DT,T)
    );

Add the case

::

    $ cp -rf porous_011_porousSimpleCopyDarcyBCFoam porous_012_porousSimpleTemperatureFoam
    $ cd porous_012_porousSimpleTemperatureFoam
    $ rm -rf 0.* [1-9]* postProcessing *.foam 	    
    
    
Create temperature boundary conditions:

::

    $ cp 0/p 0/T
    $ nano 0/T
    
::

    FoamFile
    {
            version     2.0;
            format      ascii;
            class       volScalarField;
            object      T;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensions      [0 0 0 1 0 0 0];

    internalField   uniform 273;

    boundaryField
    {
        inlet
        {
            type            fixedValue;
            value	        uniform 573;
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
        
Add coefficients for temperature equation for transportProperties

Edit the transport properties

::

    $ gedit constant/transportProperties
    
    
::

    eps            eps [ 0 0 0 0 0 0 0 ] 0.4;

    DT             DT [ 1 1 -3 -1 0 0 0 ] 1e-02;

    rhoCps         rhoCps [ 1 -1 -2 -1 0 0 0 ] 2e4;

    rhoCpf         rhoCpf [ 1 -1 -2 -1 0 0 0 ] 5e3;    


Edit the solver properties add temperature:

::

    solvers
    {
        ...

        T
        {
            solver          PBiCG;
            preconditioner  DILU;
            tolerance       1e-06;
            relTol          0;
        }
    }

Edit the schemes:

::

	divSchemes
	{
		...
		div(phi,T) bounded Gauss linear;
	}

	laplacianSchemes
	{
		...
		laplacian(DT,T) Gauss linear corrected;
	}

	fluxRequired
	{
		...
		T               ;
	}


Add a probe to controlDict and change application.

::

	application     porousSimpleTemperatureFoam;

	functions
	{
		probes
		{
		    type            probes;
		    functionObjectLibs      ("libsampling.so");
		    enabled         true;
		    outputControl   timeStep;
		    outputInterval  1;

		    fields
		    (
		        T
		    );

		    probeLocations
		    (
		        (2 0.05 0.05) // Probe 1
		        (5 0.05 0.05) // Probe 2
		        (9 0.05 0.05) // Probe 3
		    );
		}
	}

Run the case:

::

        $ rm -rf 0.* [1-9]* postProcessing *.foam
        $ porousSimpleTemperatureFoam



Transient Solver - Simple Foam

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTemperatureFoam simpleTemperatureTransientFoam
    $ cd simpleTemperatureTransientFoam
    $ nano Make/files
        simpleTemperatureFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/simpleTemperatureFoam
    $ mv simpleTemperatureFoam.C simpleTemperatureTransientFoam.C
    $ wclean
    $ wmake

Transient Solver - Porous Simple Foam

::

    $ mv porousSimpleTemperatureFoam porousSimpleTemperatureTransientFoam
    $ cd porousSimpleTemperatureTransientFoam
    $ nano Make/files
    
        porousSimpleCopyFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/porousSimpleCopyFoam
    
    $ mv porousSimpleTemperatureFoam.C porousSimpleTemperatureTransientFoam.C
    $ wclean
    $ wmake
    
    
Add this to TEqn.H

::

    solve
    (
            (eps*rhoCpf+(1.0-eps)*rhoCps)*fvm::ddt(T) 
            + rhoCpf*fvm::div(phi,T)
            ==	
            fvm::laplacian(DT,T)
    );    
    
Remake code:
    
::
    
    $ wclean
    $ wmake   

Create Case:
    
::

    $ cp -rf porous_012_porousSimpleTemperatureFoam porous_013_porousSimpleTemperatureTransientFoam
    $ cd porous_013_porousSimpleTemperatureTransientFoam
    $ rm -rf 0.* [1-9]* postProcessing *.foam     

Edit fvSchemes for transient case   
    
::
    
    ddtSchemes
    {
        default         Euler;
    }

    divSchemes
    {
        ...
        div(phi,T) Gauss linear;
    }


Edit controlDict for transient case

::

    application     porousSimpleTemperatureTransientFoam;

    endTime         60000;

    deltaT          100;

    writeControl    runTime;

    writeInterval   1000;
    
    
Try removing the convection term in the solver:

Transient Solver - Simple Foam

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTemperatureTransientFoam simpleTemperatureTransientNoConvectionFoam
    $ cd simpleTemperatureTransientNoConvectionFoam
    $ nano Make/files
        simpleTemperatureTransientNoConvectionFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/simpleTemperatureTransientNoConvectionFoam
    $ mv simpleTemperatureTransientFoam.C simpleTemperatureTransientNoConvectionFoam.C
    $ wclean
    $ wmake

Transient Solver - Porous Simple Foam

::

    $ mv porousSimpleTemperatureTransientFoam porousSimpleTemperatureTransientNoConvectionFoam
    $ cd porousSimpleTemperatureTransientNoConvectionFoam
    $ nano Make/files
    
        porousSimpleTemperatureTransientNoConvectionFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/porousSimpleTemperatureTransientNoConvectionFoam
    
    $ mv porousSimpleTemperatureTransientFoam.C porousSimpleTemperatureTransientNoConvectionFoam.C
    $ wclean
    $ wmake

Create the case

::

    $ cp -rf porous_013_porousSimpleTemperatureTransientFoam porous_014_porousSimpleTemperatureTransientNoConvectionFoam
    $ cd porous_014_porousSimpleTemperatureTransientNoConvectionFoam
    $ rm -rf 0.* [1-9]* postProcessing *.foam

Try reducing the darcy factor and including density of one. 
Density of one = same, reducing Darcy factor = different.

::

    $ cp -rf porous_013_porousSimpleTemperatureTransientFoam porous_014_porousSimpleTemperatureTransientNoConvectionFoam
    $ cd porous_014_porousSimpleTemperatureTransientNoConvectionFoam
    $ rm -rf 0.* [1-9]* postProcessing *.foam

Try removing viscosity and convection

Transient Solver - Simple Foam

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTemperatureTransientNoConvectionFoam simpleTemperatureTransientNoConvectionDiffusionFoam
    $ cd simpleTemperatureTransientNoConvectionDiffusionFoam
    $ nano Make/files
        simpleTemperatureTransientNoConvectionDiffusionFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/simpleTemperatureTransientNoConvectionDiffusionFoam
    $ mv simpleTemperatureTransientNoConvectionFoam.C simpleTemperatureTransientNoConvectionDiffusionFoam.C
    $ wclean
    $ wmake

Transient Solver - Porous Simple Foam

::

    $ mv porousSimpleTemperatureTransientNoConvectionFoam porousSimpleTemperatureTransientNoConvectionDiffusionFoam
    $ cd simpleTemperatureTransientNoConvectionDiffusionFoam
    $ nano Make/files
    
        porousSimpleTemperatureTransientNoConvectionDiffusionFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/porousSimpleTemperatureTransientNoConvectionDiffusionFoam
    
    $ mv porousSimpleTemperatureTransientNoConvectionFoam.C porousSimpleTemperatureTransientNoConvectionDiffusionFoam.C
    $ wclean
    $ wmake

Create the case

::

    $ cp -rf porous_014_porousSimpleTemperatureTransientNoConvectionFoam porous_016_porousSimpleTemperatureTransientNoConvectionDiffusionFoam
    $ cd porous_016_porousSimpleTemperatureTransientNoConvectionDiffusionFoam
    $ rm -rf 0.* [1-9]* postProcessing *.foam

Create the case

::

    $ cp -rf porous_015_porousSimpleTemperatureTransientDarcyFoam porous_017_porousSimpleTemperatureTransientDarcyVanLeerFoam
    $ cd porous_017_porousSimpleTemperatureTransientDarcyVanLeerFoam
    $ rm -rf 0.* [1-9]* postProcessing *.foam



Change scheme for vanLeer

::

	$ nano system/fvSchemes

	divSchemes
	{
		...
		div(phi,T)		Gauss vanLeer;
	}

Run the case:

::

	$ porousSimpleTemperatureTransientFoam

Change Darcy BC so that only nu is used:
----------------------------------------


::

    $ cp -rf darcyGradPressure darcyGradPressureIncompressible
    $ cd darcyGradPressureIncompressible
    $ nano Make/files

    darcyGradPressureIncompressibleFvPatchScalarField.C

            LIB = $(FOAM_USER_LIBBIN)/ldarcyGradPressureIncompressible

    $ mv darcyGradPressureFvPatchScalarField.C darcyGradPressureIncompressibleFvPatchScalarField.C
   

::

    // We recover the value of U at the boundary:

    const fvPatchField<vector>& U = 
            patch().lookupPatchField<volVectorField, vector>(UName_);

    // Extract the dictionary from the database:
    const dictionary& transportProperties = db().lookupObject<IOdictionary> ("transportProperties");

    // Extract viscosity and permeability:
    dimensionedScalar nu(transportProperties.lookup("nu"));
    dimensionedScalar d(transportProperties.lookup("d"));

    // The pressure gradient is evaluated at the boundary with the formula:
    // n.grad(p) = -(mu/k)n.U
    // mu.value() allows the access of the value of the object mu declared as a dimensionedScalar
    // patch().nf() returns the normal vector to the patch
    gradient() = -nu.value()*d.value()*(U & patch().nf());

    // curTimeIndex_ = this->db().time().timeIndex();

    // gradient() = snGradp;
    fixedGradientFvPatchScalarField::updateCoeffs();


Test
----


::

	cp -rf porous_017_porousSimpleTemperatureTransientDarcyVanLeerFoam/ porous_018_porousSimpleTemperatureTransientDarcyVanLeerIncompressibleFoam/
	cd porous_018_porousSimpleTemperatureTransientDarcyVanLeerIncompressibleFoam
	rm -rf 0.* [1-9]* postProcessing *.foam

Change transportProperties:

::

    nu             nu [0 2 -1 0 0 0 0] 1e-05;

    d              d  [0 -2   0 0 0 0 0] 1e9;

Change controlDict

::

	libs ("ldarcyGradPressureIncompressible.so");

Run the case:

::

	$ porousSimpleTemperature

Two Temperatures
----------------

Transient Solver - Simple Foam

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTemperatureTransientFoam simpleTwoTemperaturesTransientFoam
    $ cd simpleTwoTemperaturesTransientFoam
    $ nano Make/files
        simpleTwoTemperaturesTransientFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/simpleTwoTemperaturesTransientFoam
    $ mv simpleTemperatureTransientFoam.C simpleTwoTemperaturesTransientFoam.C
    $ wclean
    $ wmake

Transient Solver - Porous Simple Foam

::

    $ mv porousSimpleTemperatureTransientFoam porousSimpleTwoTemperaturesTransientFoam
    $ cd porousSimpleTwoTemperaturesTransientFoam
    $ nano Make/files
    
        porousSimpleTwoTemperaturesTransientFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/porousSimpleTwoTemperaturesTransientFoam
    
    $ mv porousSimpleTemperatureTransientFoam.C porousSimpleTwoTemperaturesTransientFoam.C
    $ wclean
    $ wmake

Add the fields and input parameters:


::

    Info<< "Reading field Ts\n" << endl;
    volScalarField Ts
    (
        IOobject
        (
            "Ts",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field Tf\n" << endl;
    volScalarField Tf
    (
        IOobject
        (
            "Tf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

The defintion of the diffusivity:

::

    Info<< "Reading diffusivity DTf\n" << endl;

    dimensionedScalar DTf
    (
        transportProperties.lookup("DTf")
    );

    Info<< "Reading diffusivity DTs\n" << endl;

    dimensionedScalar DTs
    (
        transportProperties.lookup("DTs")
    );


Add the temperature equations in TEqn.H

::

    fvScalarMatrix TfEqn
    (
          eps*rhoCpf*fvm::ddt(Tf)
        + rhoCpf*fvm::div(phi,Tf)
        ==
          fvm::laplacian(eps*DTf,Tf,"laplacian(DT,T)")
        - fvm::Sp(h,Tf)
        + h*Ts
    );
    TfEqn.solve();

    fvScalarMatrix TsEqn
    (
          (1.0-eps)*rhoCps*fvm::ddt(Ts)
        ==
          fvm::laplacian((1.0-eps)*DTs,Ts,"laplacian(DT,T)")
        - fvm::Sp(h,Ts)
        + h*Tf
    );
    TsEqn.solve();

Test Case
---------

::

	cp -rf porous_018_porousSimpleTemperatureTransientDarcyVanLeerIncompressibleFoam porous_019_porousSimpleTwoTemperaturesTransientFoam 
	cd porous_019_porousSimpleTwoTemperaturesTransientFoam
	rm -rf 0.* [1-9]* postProcessing *.foam
	cp T Tf
	mv T Ts

Change temperature Tf:

::

	dimensions      [0 0 0 1 0 0 0];

	internalField   uniform 573;

	boundaryField
	{
		inlet
		{
		    type            fixedValue;
		    value           uniform 273;
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

Change temperature Ts:

::

	dimensions      [0 0 0 1 0 0 0];

	internalField   uniform 573;

	boundaryField
	{
		inlet
		{
		    type            zeroGradient;
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

Change transport properties:

::


    transportModel  Newtonian;

    nu             nu [0 2 -1 0 0 0 0] 1e-05;

    // For Darcy BC:

    d              d  [0 -2   0 0 0 0 0] 1e9;

    // For temperature equation:

    eps            eps [ 0 0 0 0 0 0 0 ] 0.4;

    DTf            DTf [ 1 1 -3 -1 0 0 0 ] 1e-04;

    DTs            DTs [ 1 1 -3 -1 0 0 0 ] 1e-02;

    rhoCps         rhoCps [ 1 -1 -2 -1 0 0 0 ] 2e4;

    rhoCpf         rhoCpf [ 1 -1 -2 -1 0 0 0 ] 5e3;

    h              h [ 1 -1 -3 -1 0 0 0 ] 5e-1;

Change controlDict:

::

	application     porousSimpleTwoTemperaturesTransientFoam;

	functions
	{
		probes
		{
		    type            probes;
		    functionObjectLibs      ("libsampling.so");
		    enabled         true;
		    outputControl   timeStep;
		    outputInterval  1;

		    fields
		    (
		        Tf
		        Ts
		    );

		    probeLocations
		    (
		        (2 0.05 0.05) // Probe 1
		        (5 0.05 0.05) // Probe 2
		        (9 0.05 0.05) // Probe 3
		    );
		}
	}


Edit fvSolution:

::

	{
	...
		Ts
		{
		    solver          PCG;
		    preconditioner  DIC;
		    tolerance       1e-06;
		    relTol          0;
		}

		Tf
		{
		    solver          PBiCG;
		    preconditioner  DILU;
		    tolerance       1e-06;
		    relTol          0;
		}
	}

	SIMPLE
	{
		nNonOrthogonalCorrectors 2;
	}

Change fvSchemes:

::

	divSchemes
	{
		...
		div(phi,Tf) Gauss vanLeer;
	    div(phi,Ts) Gauss vanLeer;
	}

Run the case:

::

	$ porousSimpleTwoTemperaturesTransientFoam


Problem with oscillations - include relaxation in temperature 
-------------------------------------------------------------

simpleTwoTemperaturesTransientRelaxFoam

Transient Solver - Simple Foam

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTwoTemperaturesTransientFoam simpleTwoTemperaturesTransientRelaxFoam
    $ cd simpleTwoTemperaturesTransientRelaxFoam
    $ nano Make/files
        simpleTwoTemperaturesTransientRelaxFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/simpleTwoTemperaturesTransientRelaxFoam
    $ mv simpleTwoTemperaturesTransientFoam.C simpleTwoTemperaturesTransientRelaxFoam.C
    $ wclean
    $ wmake

Transient Solver - Porous Simple Foam

::

    $ mv porousSimpleTwoTemperaturesTransientFoam porousSimpleTwoTemperaturesTransientRelaxFoam
    $ cd porousSimpleTwoTemperaturesTransientRelaxFoam
    $ nano Make/files
    
        porousSimpleTwoTemperaturesTransientRelaxFoam.C
        
        EXE = $(FOAM_USER_APPBIN)/porousSimpleTwoTemperaturesTransientRelaxFoam
    
    $ mv porousSimpleTwoTemperaturesTransientFoam.C porousSimpleTwoTemperaturesTransientRelaxFoam.C
    $ wclean
    $ wmake



Add the temperature equations in TEqn.H

::

    fvScalarMatrix TfEqn
    (
          eps*rhoCpf*fvm::ddt(Tf)
        + rhoCpf*fvm::div(phi,Tf)
        ==
          fvm::laplacian(eps*DTf,Tf,"laplacian(DT,T)")
        - fvm::Sp(h,Tf)
        + h*Ts
    );
    TfEqn.relax();
    TfEqn.solve();

    fvScalarMatrix TsEqn
    (
          (1.0-eps)*rhoCps*fvm::ddt(Ts)
        ==
          fvm::laplacian((1.0-eps)*DTs,Ts,"laplacian(DT,T)")
        - fvm::Sp(h,Ts)
        + h*Tf
    );
    TfEqn.relax();
    TsEqn.solve();

Move #createFields.H in porousSimpleTwoTemperaturesTransientRelaxFoam.C


Test Case
---------

::

	cp -rf porous_019_porousSimpleTwoTemperaturesTransientFoam  porous_020_porousSimpleTwoTemperaturesTransientRelaxFoam 
	cd porous_020_porousSimpleTwoTemperaturesTransientRelaxFoam
	rm -rf 0.* [1-9]* postProcessing *.foam

Relaxation solver:

::

	application     porousSimpleTwoTemperaturesTransientRelaxFoam;

Cause of instability is the porosity was 0.15x109 in porosityProperties. This did not match the porosity in the transportProperties. vanLeer scheme can be used.

Is it because the coefficients were different or is it because the coefficients were too small?

Using 1x108 for both works - so it must have been because the coefficients were different. Hence the probelm is with the boundary condition

Darcy-Forchheimer boundary condition
------------------------------------

Copy darcyGradPressureIncompressible boundary condition

replace darcyGradPressureIncompressible with darcyForchheimer

Change the .C file lines:

dimensionedScalar f(transportProperties.lookup("f"));

gradient() = -(nu.value()*d.value() + 0.5*f.value()*fabs(U & patch().nf()))*(U & patch().nf());


Test
----

copied old folder

deleted files

added forchheimer to porosityProperties (0.1)

changed controlDict so that new BC is applied

changed 0/p so darcyForchheimer is applied


Darcy-Brinkmann-Forchheimer boundary condition
----------------------------------------------

Copy darcyForchheimer boundary condition - change to darcyBrinkmannForchheimer

replace darcyForchheimer with darcyBrinkmannForchheimer

Change the .C file lines:

dimensionedScalar eps(transportProperties.lookup("eps"));

gradient() = (turbulence->divDevReff(U)/eps.value()) - (nu.value()*d.value() + 0.5*f.value()*(U & patch().nf()))*(U & patch().nf());

change Make/files

DIDN'T WORK - USE DARCY_FORCHHEIMER

Add walls at constant temperature
---------------------------------

::

	cp -rf porous_021_porousSimpleTwoTemperaturesTransientRelaxFoamForchheimer  porous_022_porousSimpleTwoTemperaturesWalls 
	cd porous_022_porousSimpleTwoTemperaturesWalls
	rm -rf 0.* [1-9]* postProcessing *.foam

::

	constant/polyMesh/blockMeshDict

::

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

Tf

::

   walls
    {
        type            fixedValue;
        value           uniform 573;
    }

Ts

::

   walls
    {
        type            fixedValue;
        value           uniform 573;
    }


::

	$ blockMesh

::

	$ porousSimpleTwoTemperaturesTransientRelaxFoam



Convert to 2D
-------------

::

	cp -rf porous_022_porousSimpleTwoTemperaturesWalls porous_023_porousSimpleTwoTemperaturesWalls2D
	cd porous_023_porousSimpleTwoTemperaturesWalls2D
	rm -rf 0.* [1-9]* postProcessing processor* *.foam logs

::

	Lx 10;
	Ly 1;
	Lz 0.1;


	// Mesh definition (homogeneous grid with a single cell in the y and z axis since the simulation is$

	blocks
	(
		hex (0 1 2 3 4 5 6 7)
		porosity (60 10 1) simpleGrading (1 1 1)
	);

::

	$ blockMesh

::

	$ porousSimpleTwoTemperaturesTransientRelaxFoam



Now add the constant heat flux boundary condition
-------------------------------------------------

A) change the solver to deal with single temp

::
    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTwoTemperaturesTransientRelaxFoam simpleOneTemperatureTransientRelaxFoam
    $ cd simpleOneTemperatureTransientRelaxFoam

    $ mv simpleTwoTemperaturesTransientRelaxFoam.C simpleOneTemperatureTransientRelaxFoam.C
    $ sed -i s/simpleTwoTemperaturesTransientRelaxFoam/simpleOneTemperatureTransientRelaxFoam/g Make/files
    $ wclean
    $ wmake

1) edit createFields.H:

Declare only one temp field
Read onoly one diffusivity

::

    $ gedit createFields.H

            #include "createPhi.H"

            Info<< "Reading field T\n" << endl;
            volScalarField T
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            );


            Info<< "Reading diffusivity DT\n" << endl;

            dimensionedScalar DT
            (
                transportProperties.lookup("DT")
            );

2) Edit the sub-directory

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/simpleOneTemperatureTransientRelaxFoam
    $ mv porousSimpleTwoTemperaturesTransientRelaxFoam porousSimpleOneTemperatureTransientRelaxFoam
    $ cd porousSimpleOneTemperatureTransientRelaxFoam
    $ sed -i s/porousSimpleTwoTemperaturesTransientRelaxFoam/porousSimpleOneTemperatureTransientRelaxFoam/g Make/files
    $ mv porousSimpleTwoTemperaturesTransientRelaxFoam.C porousSimpleOneTemperatureTransientRelaxFoam.C

Edit the TEqn.H file

::

    $ gedit TEqn.H


            fvScalarMatrix TEqn
            (

                        (eps*rhoCpf+(1.0-eps)*rhoCps)*fvm::ddt(T)
                + rhoCpf*fvm::div(phi,T)
                ==
                    fvm::laplacian(DT,T)
            );
            
            TEqn.relax();
            TEqn.solve();


    $ wclean
    $ wmake

3) test the new solver

::

    $ cp -rf porous_023_porousSimpleTwoTemperaturesWalls2D porous_024_porousSimpleOneTemperatureWalls2D 
    $ cd porous_024_porousSimpleOneTemperatureWalls2D
    $ rm -rf 0.* [1-9]* postProcessing *.foam
    $ cd 0
    $ rm Ts
    $ mv Tf T
    $ nano T

            object T;

    $ nano system/controlDict

            application     porousSimpleOneTemperatureTransientRelaxFoam;

    fields
    (
        T
    );

    $ nano 0/T

            internalField   uniform 273;

            inlet
            {
                type            fixedValue;
                    value		uniform 573;
            }

    $ nano constant/transportProperties

            DT             DT [ 1 1 -3 -1 0 0 0 ] 1e-02;

    $ nano system/fvSolution

            T 
            {
                solver          PBiCG;
                preconditioner  DILU;
                tolerance       1e-06;
                relTol          0;
            }

    $ nano system/fvSchemes

            divSchemes
            {
                    div(phi,T) Gauss vanLeer;
            }

            laplacianSchemes
            {

                    laplacian((k|mu),p) Gauss linear corrected;

            }

            fluxRequired
            {
                    default         no;
                    p               ;
                T               ;
            }


Run the solver - make sure it is the correct solver!!!!! not the parent solver

::

	$ porousSimpleOneTemperatureTransientRelaxFoam

Edit the gnuplot file:

::

    $ nano plot_probes

    set ylabel "temperature (K)"
    set xlabel "time (s)"

    plot "postProcessing/probes/0/T" using 1:2 with lines lw 4 linecolor rgb "red" title "T Probe x=2m" ,\
        "postProcessing/probes/0/T" using 1:3 with lines lw 4 linecolor rgb "green" title "T Probe x=5m" ,\
        "postProcessing/probes/0/T" using 1:4 with lines lw 4 linecolor rgb "blue" title "T Probe x=9m"

    $ gnuplot

4) Change the walls so that upper_wall and lower_wall is created

copy the previous case directory and rename

Edit the blockMeshDict to change walls to wall_one and wall_two

::

	$ gedit constant/polyMesh/blockMeshDict 

		wall_upper
		{
		    type wall;
		    faces
		    (
		        (3 7 6 2)
		    );
		}

		wall_lower
		{
		    type wall;
		    faces
		    (
		        (0 1 5 4)
		    );
		}

	$ blockMesh


Change the boundary conditions

::

	$ gedit 0/p

		wall_upper
		{
			type		darcyForchheimer;
			value		uniform 0;
		}

		wall_lower
		{
			type		darcyForchheimer;
			value		uniform 0;
		}

	$ gedit 0/T

		wall_upper
		{
			type		fixedValue;
			value		uniform 273;
		}

		wall_lower
		{
			type		fixedValue;
			value		uniform 273;
		}

	$ gedit 0/U

		wall_upper
		{
			type		fixedValue;
			value		uniform (0 0 0);
		}

		wall_lower
		{
			type		fixedValue;
			value		uniform (0 0 0);
		}

Run the solver - make sure it is the correct solver!!!!! not the parent solver

::

	$ porousSimpleOneTemperatureTransientRelaxFoam


B) add a constant heat flux boundary condition

::

	$ cp -rf porous_025_porousSimpleOneTemperatureWalls2DUpperLower porous_026_porousSimpleOneTemperatureWalls2DUpperLowerHeatFlux   
	$ cd porous_026_porousSimpleOneTemperatureWalls2DUpperLowerHeatFlux
	$ rm -rf 0.* [1-9]* postProcessing *.foam *.png

make the boundary condition:

::

	$ cd $WM_PROJECT_USER_DIR/boundary/wallHeatFluxIncompressible
	$ wclean libso
	$ wmake libso



add the library:

::

	$ nano system/controlDict

		libs 
		(
			"ldarcyForchheimer.so"
			"libwallHeatFluxIncompressible.so"
		);


Run the solver - make sure it is the correct solver!!!!! not the parent solver

	$ porousSimpleOneTemperatureTransientRelaxFoam

C) add the constant heat flux condition to the lower wall


::

	$ nano 0/T

		wall_lower
		{
			type		wallHeatFluxIncompressible;
			value		uniform 273;
		}


Parabolic velocity field
------------------------

::

	$ mkdir -p $WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/derived/parabolicVelocity
	$ cd $WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/derived/parabolicVelocity

Download .C and .H files from: https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-Core-OpenFOAM-1.5-dev/tree/master/src/finiteVolume/fields/fvPatchFields/derived/parabolicVelocity

    1. Go to the file you want to download.
    2. Click it to view the contents within the GitHub UI.
    3. In the top right, right click the Raw button.
    4. Save as...

::

	$ mkdir Make
	$ cd Make
	$ gedit files

		parabolicVelocityFvPatchVectorField.C
		LIB = $(FOAM_USER_LIBBIN)/libmyBCs

	$ gedit options

		EXE_INC = \
		-I$(LIB_SRC)/finiteVolume/lnInclude
		EXE_LIBS =

::

	$ run
	$ cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily porous_029_pitzDailyParabolicInletDynamicLibrary
	$ cd porous_029_pitzDailyParabolicInletDynamicLibrary
	$ blockMesh

::

    $ nano 0/U

            inlet
            {
                    type            parabolicVelocity;
                    n               (1 0 0);
                    y               (0 1 0);
                    maxValue        1;
                    value           uniform (0 0 0); // Dummy for paraFoam
            }

::

    $ system/controlDict

            libs ("libmyBCs.so");

Try this on the porous medium case
----------------------------------

::

    $ nano 0/U

            inlet
            {
                    type            parabolicVelocity;
                    n               (1 0 0);
                    y               (0 1 0);
                    maxValue        1;
                    value           uniform (0 0 0); // Dummy for paraFoam
            }

::

    $ system/controlDict

            libs ("libmyBCs.so");


Alter velocity profile by changing the formula, without introducing new variables
---------------------------------------------------------------------------------

Create a new folder constantCoefficientsProfile


::
	
	$ wclean libso
	$ mv parabolicVelocityFvPatchVectorField.C constantCoefficientsFvPatchVectorField.C
	$ mv parabolicVelocityFvPatchVectorField.H constantCoefficientsFvPatchVectorField.H
	$ sed -i s/parabolicVelocity/constantCoefficients/g constantCoefficientsFvPatchVectorField.C
	$ sed -i s/parabolicVelocity/constantCoefficients/g constantCoefficientsFvPatchVectorField.H
	$ sed -i s/parabolicVelocity/constantCoefficients/g Make/files


Test the solver, changing:


::

	$ system/controlDict

		libs 
		(
			"ldarcyForchheimer.so"
			"libconstantCoefficients.so"
		);

	$ 0/U

		inlet
		{
		    type            constantCoefficients;
		    n               (1 0 0);
		    y               (0 1 0);
		    maxValue        1;
		    value           uniform (0 0 0); // Dummy for paraFoam
		}



Alter velocity profile by changing the formula, with introducing new variables
------------------------------------------------------------------------------

Create a new folder constantCoefficientsProfile


::
	
	$ wclean libso
	$ mv constantCoefficientsFvPatchVectorField.C hyperbolicVelocityFvPatchVectorField.C
	$ mv constantCoefficientsFvPatchVectorField.H hyperbolicVelocityFvPatchVectorField.H
	$ sed -i s/constantCoefficients/hyperbolicVelocity/g hyperbolicVelocityFvPatchVectorField.C
	$ sed -i s/constantCoefficients/hyperbolicVelocity/g hyperbolicVelocityFvPatchVectorField.H
	$ sed -i s/constantCoefficients/hyperbolicVelocity/g Make/files


Test the solver, changing:


::

	$ system/controlDict

		libs 
		(
			"ldarcyForchheimer.so"
			"libhyperbolicVelocity.so"
		);

	$ 0/U

		inlet
		{
		    type            hyperbolicVelocity;
		    n               (1 0 0);
		    y               (0 1 0);
		    maxValue        1;
		    value           uniform (0 0 0); // Dummy for paraFoam
		}



Introduce Brinkmann term
------------------------

::

    $ cp -rf porous_034_porousSimpleTwoTemperaturesWalls2DPorousHyperbolicInletDarcyForchheimerCalculated porous_035_porousSimpleTwoTemperaturesWalls2DHyperbolicBrinkmann   
    $ cd porous_035_porousSimpleTwoTemperaturesWalls2DHyperbolicBrinkmann


::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume
    $ mkdir cfdTools
    $ cd cfdTools
    $ mkdir general
    $ cd general
    $ mkdir porosityModel
    $ mkdir $WM_PROJECT_USER_DIR/src/finiteVolume/cfdTools/general/porosityModel/
    $ mkdir /cfdTools/general/porosityModel/
    $ cp -r $FOAM_SRC/finiteVolume/cfdTools/general/porosityModel/DarcyForchheimer DarcyBrinkmannForchheimer
    
    
::
    
    $ cd DarcyBrinkmannForchheimer
    $ mv DarcyForchheimer.C DarcyBrinkmannForchheimer.C
    $ mv DarcyForchheimer.H DarcyBrinkmannForchheimer.H
    $ mv DarcyForchheimerTemplates.C DarcyBrinkmannForchheimerTemplates.C
    
::

    $ sed -i s/DarcyForchheimer/DarcyBrinkmannForchheimer/g DarcyBrinkmannForchheimer.C
    $ sed -i s/DarcyForchheimer/DarcyBrinkmannForchheimer/g DarcyBrinkmannForchheimer.H
    $ sed -i s/DarcyForchheimer/DarcyBrinkmannForchheimer/g DarcyBrinkmannForchheimerTemplates.C
     
    
::

    $ cd $WM_PROJECT_USER_DIR/finiteVolume
    $ mkdir Make
    $ cd Make
    $ gedit files

    general = cfdTools/general
    porosity = $(general)/porosityModel
    $(porosity)/DarcyBrinkmannForchheimer/DarcyBrinkmannForchheimer.C
    LIB = $(FOAM_USER_LIBBIN)/libmyfiniteVolume

    $ gedit options

            EXE_INC = \
    -I$(LIB_SRC)/triSurface/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/finiteVolume/lnInclude
    LIB_LIBS = \
    -lOpenFOAM \
    -ltriSurface \
    -lmeshTools
    
::

    $ wclean libso
    $ wmake libso


Add library to system/controlDict:

::

    $ system/controlDict

            libs 
            (
                    "ldarcyForchheimer.so"
                    "libhyperbolicVelocity.so"
                    "libmyfiniteVolume.so"
            );

::

    $ constant/porosityProperties

        type            DarcyBrinkmannForchheimer;
        active          yes;
        cellZone        porosity;

        DarcyBrinkmannForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (1.0e9 0 0);
            f   f [0 -1 0 0 0 0 0] (0.1 0 0);

            coordinateSystem
            {
                type    cartesian;
                origin  (0 0 0);
                coordinateRotation
                {
                    type    axesRotation;
                    e1      (1 0 0);
                    e2      (0 1 0);
                }
            }
        }



Remove transient component - porousSimpleTwoTemperaturesSteadyRelaxFoam
-----------------------------------------------------------------------

::
    $ mv simpleTwoTemperaturesTransientRelaxFoam.C simpleTwoTemperaturesSteadyRelaxFoam.C
    $ mv porousSimpleTwoTemperaturesTransientRelaxFoam.C porousSimpleTwoTemperaturesSteadyRelaxFoam.C


::

    fvScalarMatrix TfEqn
        (
              rhoCpf*fvm::div(phi,Tf)
            ==
              fvm::laplacian(eps*DTf,Tf,"laplacian(DT,T)")
            - fvm::Sp(h,Tf)
            + h*Ts
        );
        TfEqn.relax();
        TfEqn.solve();

        fvScalarMatrix TsEqn
        (
              fvm::laplacian((1.0-eps)*DTs,Ts,"laplacian(DT,T)")
            - fvm::Sp(h,Ts)
            + h*Tf
        );
        TsEqn.relax();
        TsEqn.solve();

Make/files

::

    simpleTwoTemperaturesSteadyRelaxFoam.C

    EXE = $(FOAM_USER_APPBIN)/simpleTwoTemperaturesSteadyRelaxFoam

::

    $ cd simpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake

Make/files

::

    porousSimpleTwoTemperaturesSteadyRelaxFoam.C

    EXE = $(FOAM_USER_APPBIN)/porousSimpleTwoTemperaturesSteadyRelaxFoam


::

    $ cd porousSimpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake


Test solver
-----------

::

    $ system/fvSchemes

    ddt {

    default none

    }

Edit system/controlDict:

::

    simpleTwoTemperaturesSteadyRelaxFoam

Edit system/sampleDict:

::

    sets
    (
        lineX1
        {
            type        midPoint;
            axis        distance;

            //- cavity. Slightly perturbed so not to align with face or edge.
            start       (0.000 0.010 0.000);
            end         (0.100 0.010 0.000);
        }

    );


Run sample utility:

::

    $ sample


::

    $ gnuplot -persist plot_temperature
    $ gnuplot -persist plot_pressure
    $ gnuplot -persist plot_temperature



Read nu from transportProperties

::

    Info<< "Reading kinematic viscosity nu\n" << endl;

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );


recompile:

::

    $ cd simpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake

::

    $ cd porousSimpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake


Replace turbulence with laminar in porousSimpleTwoTemperaturesSteadyRelaxFoam/UEqn.H and simpleTwoTemperaturesSteadyRelaxFoam/UEqn.H

::

    // Construct the Momentum equation
    Info<< "Construct the Momentum equation = " << endl;
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(nu,U)
      ==
        fvOptions(U)
    );


recompile:

::

    $ cd simpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake

::

    $ cd porousSimpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake


Divide nu by eps in porousSimpleTwoTemperaturesSteadyRelaxFoam/UEqn.H and simpleTwoTemperaturesSteadyRelaxFoam/UEqn.H

::

    // Construct the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(nu/eps,U)
      ==
        fvOptions(U)
    );


recompile:

::

    $ cd simpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake

::

    $ cd porousSimpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake


Correct Temperature Equation
----------------------------

Read in surface area density and correct notations

::

    Info<< "Reading effective fluid thermal conductivity kfe\n" << endl;

    dimensionedScalar kfe
    (
        transportProperties.lookup("kfe")
    );

    Info<< "Reading effective solid thermal conductivity kse\n" << endl;

    dimensionedScalar kse
    (
        transportProperties.lookup("kse")
    );

    Info<< "Reading fluid density times heat capacity rhoCpf\n" << endl;

    dimensionedScalar rhoCpf
    (
        transportProperties.lookup("rhoCpf")
    );

    Info<< "Reading interfacial heat transfer coefficient hsf\n" << endl;

    dimensionedScalar hsf
    (
        transportProperties.lookup("hsf")
    );

    Info<< "Reading surface area density asf\n" << endl;

    dimensionedScalar asf
    (
        transportProperties.lookup("asf")
    );

::

    // Tf treated as implicit in TfEqn, Ts treated as explicit
    fvScalarMatrix TfEqn
    (
          rhof*cpf*fvm::div(phi,Tf)
        ==
          fvm::laplacian(kfe,Tf,"laplacian(DT,T)")
        - fvm::Sp(hsf*asf,Tf)
        + hsf*asf*Ts
    );
    TfEqn.relax();
    TfEqn.solve();

    // Ts treated an implicit in TsEqn, Tf treated as explicit
    fvScalarMatrix TsEqn
    (
          fvm::laplacian(kse,Ts,"laplacian(DT,T)")
        - fvm::Sp(hsf*asf,Ts)
        + hsf*asf*Tf
    );
    TsEqn.relax();
    TsEqn.solve();

::

    // For temperature equation:

    eps            eps [ 0 0 0 0 0 0 0 ] 0.4;

    kfe            kfe [ 1 1 -3 -1 0 0 0 ] 1e-04;

    kse            kse [ 1 1 -3 -1 0 0 0 ] 1e-04;

    rhof           rhof [ 1 -3 0 0 0 0 0 ] 1.184;

    cpf            cpf [ 0 2 -2 -1 0 0 0 ] 1007.0;

    hsf	           hsf [ 1 0 -3 -1 0 0 0 ] 5e-1;

    asf			   asf [ 0 -1 0 0 0 0 0 ] 1.0;


recompile:

::

    $ cd simpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake

::

    $ cd porousSimpleTwoTemperaturesSteadyRelaxFoam
    $ wclean
    $ wmake




Add fixedMeanValue pressure to outlet boundary.
-----------------------------------------------

::

	$ mkdir -p $WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/derived/fixedMeanValue
	$ cd $WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/derived/fixedMeanValue

Download .C and .H files from: https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-Core-OpenFOAM-1.5-dev/tree/master/src/finiteVolume/fields/fvPatchFields/derived/fixedMeanValue

    1. Go to the file you want to download.
    2. Click it to view the contents within the GitHub UI.
    3. In the top right, right click the Raw button.
    4. File> Save Page as...


::

	$ mkdir Make
	$ cd Make
	$ gedit files

		fixedMeanValueFvPatchField.C
		LIB = $(FOAM_USER_LIBBIN)/libfixedMeanValue

	$ gedit options

		EXE_INC = \
		-I$(LIB_SRC)/finiteVolume/lnInclude
		EXE_LIBS =

	$ cd ..

	$ wclean libso
	$ wmake libso


!!!!!!!!!!DIDN'T WORK!!!!!!!!!!!!!!!!!


Instead use fixedMean:

::

    outlet
    {
        type            fixedMean;
        meanValue       0.0;
        value           uniform 0;
    }


Still no influence on profile, so the mass balance or the grid must be having an influence.

Adjust profile parameters to match mean velocity profile


Alter velocity profile by changing the formula, with introducing new variables
------------------------------------------------------------------------------

Create a new folder hyperbolicVelocityCalculated

::
	$ mkdir hyperbolicVelocityCalculated
	$ cp -rf hyperbolicVelocity hyperbolicVelocityCalculated


::
	
	$ wclean libso
	$ mv hyperbolicVelocityFvPatchVectorField.C hyperbolicVelocityCalculatedFvPatchVectorField.C
	$ mv hyperbolicVelocityFvPatchVectorField.H hyperbolicVelocityCalculatedFvPatchVectorField.H
	$ sed -i s/hyperbolicVelocity/hyperbolicVelocityCalculated/g hyperbolicVelocityCalculatedFvPatchVectorField.C
	$ sed -i s/hyperbolicVelocity/hyperbolicVelocityCalculated/g hyperbolicVelocityCalculatedFvPatchVectorField.H
	$ sed -i s/hyperbolicVelocity/hyperbolicVelocityCalculated/g Make/files


1) Add private member variables to hyperbolicVelocityCalculatedFvPatchVectorField.H:

::

        //- Porosity
        scalar epsilon_;

        //- Forchheimer Coefficient
        scalar forch_;

        //- Plate height
        scalar height_;

        //- Fluid kinematic viscosity
        scalar kinVisc_;

        //- Permeability
        scalar perm_;


2) Add private member functions to hyperbolicVelocityCalculatedFvPatchVectorField.H:

::

        //- Porosity
        scalar& epsilon()
        {
            return epsilon_;
        }

        //- Forchheimer Coefficient
        scalar& forch()
        {
            return forch_;
        }

        //- height
        scalar& height()
        {
            return height_;
        }

        //- height
        scalar& kinVisc()
        {
            return kinVisc_;
        }

        //- Permeability
        scalar& perm()
        {
            return perm_;
        }

Recompile:

::

	$ wclean libso
	$ wmake libso



3) Create default values in hyperbolicVelocityCalculatedFvPatchVectorField.C:

::

	epsilon_(1),
	forch_(1),
	height_(1),
	kinVisc_(1),
	perm_(1)

::

    epsilon_(ptf.epsilon_),
    forch_(ptf.forch_),
    height_(ptf.height_,
    kinVisc_(ptf.kinVisc_),
    perm_(ptf.perm_)

Lookup the values:

::

	epsilon_(dict.lookup("epsilon")),
	forch_(dict.lookup("forch")),
	height_(dict.lookup("height")),
	kinVisc_(dict.lookup("kinVisc")),
	perm_(dict.lookup("perm"))

Add this code:

::

        const double Lambda = ((pow(epsilon_, (3.0/2.0)))*forch_*maxValue_*height_) / kinVisc_ ;
        const double Darcy = perm_ / (sqr(height_)*epsilon_);
        const double ACoeff = (2.0 / 3.0)*Lambda*(pow(Darcy, (-1.0/2.0)));
        const double BCoeff = (pow(Darcy, (-1.0))) + (4.0/3.0)*Lambda*pow(Darcy, (-1.0/2.0));
        const double DCoeff = (sqrt(ACoeff + BCoeff))/2.0; 
        const double Value = sqrt(ACoeff/(ACoeff + BCoeff));
        const double CCoeff = -(1.0 / DCoeff)*acosh(1.0/Value) - 1.0;

        vectorField::operator=(n_*maxValue_*(1.0 - ((ACoeff+BCoeff)/ACoeff)*sqr(  1.0 / (cosh (DCoeff*  (mag(coord) + CCoeff) ) )   )));

        Info << "Lambda" << Lambda << endl;
        Info << "Darcy" << Darcy << endl;
        Info << "ACoeff" << ACoeff << endl;
        Info << "BCoeff" << BCoeff << endl;
        Info << "DCoeff" << DCoeff << endl;
        Info << "CCoeff" << CCoeff << endl;


Recompile:

::

	$ wclean libso
	$ wmake libso

**changed so that Darcy and Forchheimer coefficients are calculated** 

Test the solver, changing:

::

	$ system/controlDict

		libs 
		(
			"ldarcyForchheimer.so"
			"libhyperbolicVelocityCalculated.so"
		);

	$ 0/U

		inlet
		{
	    type            hyperbolicVelocityCalculated;
	    n               (1 0 0);
	    y               (0 1 0);
	    maxValue        1.0e-4;
		epsilon			0.380;
		diameter		0.002;
		height			0.010;
		kinVisc 		1.0e-5;
	    value           uniform (0 0 0); // Dummy for paraFoam
		}




**changed so that Darcy and Forchheimer coefficients are calculated** 

Checked coefficients against Python code - coefficients were the same

Now make the DarcyBrinkmannForchheimer into DarcyForcheimerCalculated

::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/
    $ wclean libso
    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/cfdTools/general/porosityModel

    $ cp -rf DarcyBrinkmannForchheimer DarcyForchheimerCalculated
    $ cd DarcyForchheimerCalculated

    $ mv DarcyBrinkmannForchheimer.C DarcyForchheimerCalculated.C
    $ mv DarcyBrinkmannForchheimer.H DarcyForchheimerCalculated.H
    $ mv DarcyBrinkmannForchheimerTemplates.C DarcyForchheimerCalculatedTemplates.C
    
::

    $ sed -i s/DarcyBrinkmannForchheimer/DarcyForchheimerCalculated/g DarcyForchheimerCalculated.C
    $ sed -i s/DarcyBrinkmannForchheimer/DarcyForchheimerCalculated/g DarcyForchheimerCalculated.H
    $ sed -i s/DarcyBrinkmannForchheimer/DarcyForchheimerCalculated/g DarcyForchheimerCalculatedTemplates.C



::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/Make
    $ gedit files

    general = cfdTools/general
    porosity = $(general)/porosityModel
    $(porosity)/DarcyForchheimerCalculated/DarcyForchheimerCalculated.C
    LIB = $(FOAM_USER_LIBBIN)/libmyfiniteVolume

    $ gedit options

            EXE_INC = \
    -I$(LIB_SRC)/triSurface/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/finiteVolume/lnInclude
    LIB_LIBS = \
    -lOpenFOAM \
    -ltriSurface \
    -lmeshTools


::

    $ cd ..
    $ wclean libso
    $ wmake libso




Add library to system/controlDict:

::

    $ system/controlDict

            libs 
            (
                    "ldarcyForchheimer.so"
                    "libhyperbolicVelocity.so"
                    "libmyfiniteVolume.so"
            );

::

    $ constant/porosityProperties

    type            DarcyForchheimerCalculated;
    active          yes;
    cellZone        porosity;

    DarcyForchheimerCalculatedCoeffs
    {
        d   d [0 -2 0 0 0 0 0] (1.0e9 0 0);
        f   f [0 -1 0 0 0 0 0] (0.1 0 0);

        coordinateSystem
        {
            type    cartesian;
            origin  (0 0 0);
            coordinateRotation
            {
                type    axesRotation;
                e1      (1 0 0);
                e2      (0 1 0);
            }
        }
    }



Installed SWAK4FOAM

tested using cavity case:


::

    $ nano 0/U

            movingWall
            {
                type groovyBC;
                valueExpression "vector(sin(w*time()), 0, 0)";
                variables "w=pi;";
                value $internalField;
            }

    $ nano system/controlDict

            endTime         2.0;
            deltaT          0.001;

            libs ("libgroovyBC.so");

    $ rm -rf 0.* [1-9]* postProcessing *.foam

    $ icoFoam > log


Now apply Ts = Tf using Swak4Foam:



::

    $ rm -rf 0.* [1-9]* postProcessing *.foam

    $ nano 0/Ts

            wall_upper
            {
                    type groovyBC;
                    valueExpression "Tf";
                    variables "Tf@wall_upper=Tf;";
                    value $internalField;
            } 

    $ nano system/controlDict

            libs ("libgroovyBC.so");


    $ porousSimpleTwoTemperaturesSteadyRelaxFoam > log &


There is a difference of around 20K between the surface temperature Ts and Tf.

Change the inlet so that there is constant gradient for Ts:

::

    $ nano 0/Ts

    inlet
    {
            type    zeroGradient;
    }


It could be that the temperature has not reached steady state, so the temperature applied is from the previous timestep previous time-step or iteration are utilised

"When evaluating the expression in run-time, the state of the fields on the remote patches from the previous time-step or iteration are utilised" - false, the values are not changing with time!!!

Influence of mesh is also possible


Correct Darcy-Forchheimer pressure boundary so that coefficients are calulated
------------------------------------------------------------------------------

Easier

::

    $ cp -rf darcyForchheimer darcyForchheimerCalculated
    $ cd darcyForcheimerCalculated
    $ nano Make/files
    
        darcyForchheimerCalculatedFvPatchScalarField.C

            LIB = $(FOAM_USER_LIBBIN)/ldarcyForchheimerCalculated

    $ mv darcyForchheimerFvPatchScalarField.C darcyForchheimerCalculatedFvPatchScalarField.C
    $ mv darcyForchheimerFvPatchScalarField.H darcyForchheimerCalculatedFvPatchScalarField.H


::

    $ sed -i s/darcyForchheimer/darcyForchheimerCalculated/g darcyForchheimerCalculatedFvPatchScalarField.C
    $ sed -i s/darcyForchheimer/darcyForchheimerCalculated/g darcyForchheimerCalculatedFvPatchScalarField.H

    
::

    // Extract viscosity, diameter and epsilon:
    dimensionedScalar nu(transportProperties.lookup("nu"));
    dimensionedScalar diameter(transportProperties.lookup("diameter"));
    dimensionedScalar eps(transportProperties.lookup("eps"));


    dimensionedScalar darcy =
    dimensionedScalar
    (
            "darcy", /// A name field
            dimensionSet( 0, 2, 0, 0, 0, 0, 0 ),
            1.0 /// Value to initialize
    );


    dimensionedScalar forchheimer =
    dimensionedScalar
    (
            "forchheimer", /// A name field
            dimensionSet( 0, 0, 0, 0, 0, 0, 0 ),
            1.0 /// Value to initialize
    );

    darcy = ( (pow(eps, 3.0))*(sqr(diameter)) ) / (150.0*sqr(1.0 - eps));

    forchheimer = 1.75 / ((sqrt(150.0))* pow(eps, (3.0/2.0)));

    Info << "Darcy Coefficient" << darcy.value() << endl;
    Info << "Forchheimer Coefficient" << forchheimer.value() << endl;
    
    gradient() = -(nu.value()/darcy.value() + ( forchheimer.value()/sqrt(darcy.value()) )*(U & patch().nf()))*(U & patch().nf());
    
    
::

    $ wclean libso
    $ wmake libso



The function .value() returns a dimensionless number

To have a number with dimensions use the dimensioned scalar variable


Correct Darcy-Forchheimer porosity model so that coefficients are calculated
----------------------------------------------------------------------------

1) Read in epsilon and diameter

The apply function creates the source term:

::

    template<class RhoFieldType>
    void Foam::porosityModels::DarcyForchheimer::apply
    (
        scalarField& Udiag,
        vectorField& Usource,
        const scalarField& V,
        const RhoFieldType& rho,
        const scalarField& mu,
        const vectorField& U
    ) const
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const tensorField& dZones = D_[zoneI];
            const tensorField& fZones = F_[zoneI];

            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label cellI = cells[i];
                const label j = this->fieldIndex(i);
                const tensor Cd =
                    mu[cellI]*dZones[j] + (rho[cellI]*mag(U[cellI]))*fZones[j];

                const scalar isoCd = tr(Cd);

                Udiag[cellI] += V[cellI]*isoCd;
                Usource[cellI] -= V[cellI]*((Cd - I*isoCd) & U[cellI]);
            }
        }
    }


Use Momentum Equation of Zhao
-----------------------------

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTwoTemperaturesSteadyRelaxFoam simpleTwoTemperaturesZhaoFoam
    $ cd simpleTwoTemperaturesZhaoFoam

    $ mv simpleTwoTemperaturesSteadyRelaxFoam.C simpleTwoTemperaturesZhaoFoam.C
    $ sed -i s/simpleTwoTemperaturesSteadyRelaxFoam/simpleTwoTemperaturesZhaoFoam/g Make/files

    $ nano UEqn.H

            tmp<fvVectorMatrix> UEqn
            (
                (1.0/eps)*fvm::div(phi, U)
                - fvm::laplacian(nu/eps, U)
                ==
                fvOptions(U)
            );

    $ wclean
    $ wmake

2) Edit the sub-directory

::

    $ mv porousSimpleTwoTemperaturesSteadyRelaxFoam porousSimpleTwoTemperaturesZhaoFoam
    $ cd porousSimpleOneTemperatureTransientRelaxFoam
    $ sed -i s/porousSimpleTwoTemperaturesSteadyRelaxFoam/porousSimpleTwoTemperaturesZhaoFoam/g Make/files
    $ mv porousSimpleTwoTemperaturesSteadyRelaxFoam.C porousSimpleTwoTemperaturesZhaoFoam.C

    $ nano UEqn.H

            tmp<fvVectorMatrix> UEqn
            (
                (1.0/eps)*fvm::div(phi, U)
                - fvm::laplacian(nu/eps, U)
                ==
                fvOptions(U)
            );

    $ wclean
    $ wmake


Inertia free velocity profile (Profile of Zhao)
-----------------------------------------------

::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/derived
    $ cp -rf hyperbolicVelocityCalculated hyperbolicVelocityZhao


::

	$ cd hyperbolicVelocityZhao
	$ wclean libso
	$ mv hyperbolicVelocityCalculatedFvPatchVectorField.C hyperbolicVelocityZhaoFvPatchVectorField.C
	$ mv hyperbolicVelocityCalculatedFvPatchVectorField.H hyperbolicVelocityZhaoFvPatchVectorField.H
	$ sed -i s/hyperbolicVelocityCalculated/hyperbolicVelocityZhao/g hyperbolicVelocityZhaoFvPatchVectorField.C
	$ sed -i s/hyperbolicVelocityCalculated/hyperbolicVelocityZhao/g hyperbolicVelocityZhaoFvPatchVectorField.H
	$ sed -i s/hyperbolicVelocityCalculated/hyperbolicVelocityZhao/g Make/files


::

    $ nano hyperbolicVelocityZhaoFvPatchVectorField.H

    //- Plate height
    scalar h_;

    //- Permeability
    scalar k_;

    //- height
    scalar& h()
    {
        return h_;
    }

    //- height
    scalar& k()
    {
        return k_;
    }

::

	$ nano hyperbolicVelocityZhaoFvPatchVectorField.C

		hyperbolicVelocityZhaoFvPatchVectorField::hyperbolicVelocityZhaoFvPatchVectorField
		(
			const fvPatch& p,
			const DimensionedField<vector, volMesh>& iF
		)
		:
			fixedValueFvPatchVectorField(p, iF),
			maxValue_(0),
			n_(1, 0, 0),
			y_(0, 1, 0),
			h_(1),
			k_(1)
		{}


		hyperbolicVelocityZhaoFvPatchVectorField::hyperbolicVelocityZhaoFvPatchVectorField
		(
			const hyperbolicVelocityZhaoFvPatchVectorField& ptf,
			const fvPatch& p,
			const DimensionedField<vector, volMesh>& iF,
			const fvPatchFieldMapper& mapper
		)
		:
			fixedValueFvPatchVectorField(ptf, p, iF, mapper),
			maxValue_(ptf.maxValue_),
			n_(ptf.n_),
			y_(ptf.y_),
			h_(ptf.h_),
			k_(ptf.k_)
		{}


		hyperbolicVelocityZhaoFvPatchVectorField::hyperbolicVelocityZhaoFvPatchVectorField
		(
			const fvPatch& p,
			const DimensionedField<vector, volMesh>& iF,
			const dictionary& dict
		)
		:
			fixedValueFvPatchVectorField(p, iF),
			maxValue_(readScalar(dict.lookup("maxValue"))),
			n_(dict.lookup("n")),
			y_(dict.lookup("y")),
			h_(readScalar(dict.lookup("h"))),
			k_(readScalar(dict.lookup("k")))
		{
			if (mag(n_) < SMALL || mag(y_) < SMALL)
			{
				FatalErrorIn("hyperbolicVelocityZhaoFvPatchVectorField(dict)")
				    << "n or y given with zero size not correct"
				    << abort(FatalError);
			}

			n_ /= mag(n_);
			y_ /= mag(y_);

			evaluate();
		}


		hyperbolicVelocityZhaoFvPatchVectorField::hyperbolicVelocityZhaoFvPatchVectorField
		(
			const hyperbolicVelocityZhaoFvPatchVectorField& fcvpvf,
			const DimensionedField<vector, volMesh>& iF
		)
		:
			fixedValueFvPatchVectorField(fcvpvf, iF),
			maxValue_(fcvpvf.maxValue_),
			n_(fcvpvf.n_),
			y_(fcvpvf.y_),
			h_(fcvpvf.h_),
			k_(fcvpvf.k_)
		{}


		// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

		void hyperbolicVelocityZhaoFvPatchVectorField::updateCoeffs()
		{
			if (updated())
			{
				return;
			}

			// Get range and orientation
			boundBox bb(patch().patch().localPoints(), true);

			vector ctr = 0.5*(bb.max() + bb.min());

			const vectorField& c = patch().Cf();

			// Calculate local 1-D coordinate for the parabolic profile
			scalarField coord = 2*((c - ctr) & y_)/((bb.max() - bb.min()) & y_);

			const double Darcy = sqrt(k_) / h_;
			const double Coefficient = 
		   
			vectorField::operator=(n_*maxValue_*(1.0 - ( cosh((Darcy**(-0.5))*mag(coord)) )/cosh((Darcy**(-0.5)) ) ) );


		}


		// Write
		void hyperbolicVelocityZhaoFvPatchVectorField::write(Ostream& os) const
		{
			fvPatchVectorField::write(os);
			os.writeKeyword("maxValue")
				<< maxValue_ << token::END_STATEMENT << nl;
			os.writeKeyword("n")
				<< n_ << token::END_STATEMENT << nl;
			os.writeKeyword("y")
				<< y_ << token::END_STATEMENT << nl;
			os.writeKeyword("height")
				<< h_ << token::END_STATEMENT << nl;
			os.writeKeyword("permeability")
				<< k_ << token::END_STATEMENT << nl;
			writeEntry("value", os);
		}


Recompile:

::

	$ wclean libso
	$ wmake libso



Specify K and F in porous medium
--------------------------------


::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/cfdTools/general/porosityModel
    $ cp -rf DarcyForchheimerCalculated DarcyForchheimerCoefficients 
    
    
::
    
    $ cd DarcyForchheimerCoefficients
    $ mv DarcyForchheimerCalculated.C DarcyForchheimerCoefficients.C
    $ mv DarcyForchheimerCalculated.H DarcyForchheimerCoefficients.H
    $ mv DarcyForchheimerCalculatedTemplates.C DarcyForchheimerCoefficientsTemplates.C
    
::

    $ sed -i s/DarcyForchheimerCalculated/DarcyForchheimerCoefficients/g DarcyForchheimerCoefficients.C
    $ sed -i s/DarcyForchheimerCalculated/DarcyForchheimerCoefficients/g DarcyForchheimerCoefficients.H
    $ sed -i s/DarcyForchheimerCalculated/DarcyForchheimerCoefficients/g DarcyForchheimerCoefficientsTemplates.C
     
    
::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume
    $ cd Make
    $ gedit files

    general = cfdTools/general
    porosity = $(general)/porosityModel
    $(porosity)/DarcyForchheimerCoefficients/DarcyForchheimerCoefficients.C
    LIB = $(FOAM_USER_LIBBIN)/libmyfiniteVolume
   
::

    $ cd ..
    $ wclean libso
    $ wmake libso


::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/cfdTools/general/porosityModel/DarcyForchheimerCoefficients


    $ kate DarcyForchheimerCoefficients.H

    //- Permeability        
    dimensionedScalar k_;
    
    //- Forchheimer       
    dimensionedScalar f_;



    $ kate DarcyForchheimerCoefficients.C

            if (coordSys_.R().uniform())
            {
                forAll (cellZoneIDs_, zoneI)
                {
                    D_[zoneI].setSize(1);
                    F_[zoneI].setSize(1);

                    D_[zoneI][0] = tensor::zero;
                    D_[zoneI][0].xx() = 1.0 / k_.value();
                    D_[zoneI][0].yy() = 1.0 / k_.value();
                    D_[zoneI][0].zz() = 1.0 / k_.value();

                    D_[zoneI][0] = coordSys_.R().transformTensor(D_[zoneI][0]);

                    F_[zoneI][0] = tensor::zero;
                    F_[zoneI][0].xx() = f_.value()/sqrt(k_.value());
                    F_[zoneI][0].yy() = f_.value()/sqrt(k_.value());
                    F_[zoneI][0].zz() = f_.value()/sqrt(k_.value());

                    F_[zoneI][0] = coordSys_.R().transformTensor(F_[zoneI][0]);
                }
            }
            else
            {
                forAll(cellZoneIDs_, zoneI)
                {
                    const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

                    D_[zoneI].setSize(cells.size());
                    F_[zoneI].setSize(cells.size());

                    forAll(cells, i)
                    {
                        D_[zoneI][i] = tensor::zero;
                        D_[zoneI][i].xx() = 1.0 / k_.value();
                        D_[zoneI][i].yy() = 1.0 / k_.value();
                        D_[zoneI][i].zz() = 1.0 / k_.value();

                        F_[zoneI][i] = tensor::zero;
                        F_[zoneI][i].xx() = f_.value()/sqrt(k_.value());
                        F_[zoneI][i].yy() = f_.value()/sqrt(k_.value());
                        F_[zoneI][i].zz() = f_.value()/sqrt(k_.value());
                    }

                    const coordinateRotation& R = coordSys_.R(mesh_, cells);

                    D_[zoneI] = R.transformTensor(D_[zoneI], cells);
                    F_[zoneI] = R.transformTensor(F_[zoneI], cells);
                }
            }

::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume
    $ wclean libso
    $ wmake libso



Specify K and F in boundary
---------------------------

::

    $ cd $WM_PROJECT_USER_DIR/boundary
    $ cp -rf darcyForchheimerCalculated darcyForchheimerCoefficients
    $ cd darcyForchheimerCoefficients

    $ mv darcyForchheimerCalculatedFvPatchScalarField.C darcyForchheimerCoefficientsFvPatchScalarField.C
    $ mv darcyForchheimerCalculatedFvPatchScalarField.H darcyForchheimerCoefficientsFvPatchScalarField.H


::

    $ sed -i s/darcyForchheimerCalculated/darcyForchheimerCoefficients/g darcyForchheimerCoefficientsFvPatchScalarField.C
    $ sed -i s/darcyForchheimerCalculated/darcyForchheimerCoefficients/g darcyForchheimerCoefficientsFvPatchScalarField.H
    $ sed -i s/darcyForchheimerCalculated/darcyForchheimerCoefficients/g Make/files


Recompile:

::

	$ wclean libso
	$ wmake libso


::

	$ kate darcyForchheimerCoefficientsFvPatchScalarField.C

		dimensionedScalar nu(transportProperties.lookup("nu"));
		dimensionedScalar k(transportProperties.lookup("k"));
		dimensionedScalar f(transportProperties.lookup("f"));


		// The pressure gradient is evaluated at the boundary with the formula:
		// n.grad(p) = -(mu/k)n.U
		// mu.value() allows the access of the value of the object mu declared as a dimensionedScalar
		// patch().nf() returns the normal vector to the patch
		gradient() = -(nu.value()/k.value() + ( f.value()/sqrt(k.value()) )*(U & patch().nf()))*(U & patch().nf());


Test Zhao case
--------------

Create grid:

::

	$ cd porous_037_porousSimpleTwoTemperaturesZhaoProfile
	$ rm -rf 0.* [1-9]* postProcessing *.foam

::

	// The vertices can be defined using variables Lx, Ly, Lz. It saves time to modify the size of the domain

	Lx 0.114;     // add the ends
	Ly 0.045;
	Lz 0.045;
	Dx1 0.0057;
	Dx2 0.1197;
	Dx3 0.1254;


	vertices
	(
		(0 0 0)
		($Dx1 0 0)
		($Dx1 $Ly 0)
		(0 $Ly 0)
		(0 0 $Lz)
		($Dx1 0 $Lz)
		($Dx1 $Ly $Lz)
		(0 $Ly $Lz)
		($Dx2 0 0)
		($Dx3 0 0)
		($Dx3 $Ly 0)
		($Dx2 $Ly 0)
		($Dx2 0 $Lz)
		($Dx3 0 $Lz)
		($Dx3 $Ly $Lz)
		($Dx2 $Ly $Lz)
	);

	// Mesh definition (homogeneous grid with a single cell in the y and z axis since the simulation is 1D)

	blocks
	(
		hex (0 1 2 3 4 5 6 7) 
		porosity (6 46 1) simpleGrading (1 1 1)
		hex (1 8 11 2 5 12 15 6) 
		porosity (127 46 1) simpleGrading (1 1 1)
		hex (8 9 10 11 12 13 14 15) 
		porosity (6 46 1) simpleGrading (1 1 1)
	);

	edges
	(
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
				(10 14 13 9)
		    );
		}

		// Faces orthogonal to y and z axis are defined as «empty» to specify that the simulation is 1D

		wall_upper
		{
		    type wall;
		    faces
		    (
		        (3 7 6 2)
				(2 6 15 11)
				(11 15 14 10)
		    );
		}

		wall_lower
		{
		    type wall;
		    faces
		    (
		        (0 1 5 4)
				(8 9 13 12)
		    );
		}

		wall_lower_heated
		{
		    type wall;
		    faces
		    (
				(1 8 12 5)
		    );
		}

		frontAndBack
		{
		    type empty;
		    faces
		    (
		        (0 3 2 1)
		        (4 5 6 7)
				(5 12 15 6)
				(12 13 14 15)
				(1 2 11 8)
				(8 11 10 9)
		    );
		}

	);



Change boundary condition so lower wall is heated:

P:

::

    inlet
    {
        type            darcyForchheimerCoefficients;
        value           uniform 0;
    }

    outlet
    {
        type            fixedMean;
        meanValue       0.0;
        value           uniform 0;
    }

    wall_upper
    {
        type            zeroGradient;
    }

    wall_lower_heated
    {
        type            zeroGradient;
    }

    wall_lower
    {
        type            zeroGradient;
    }

    frontAndBack
    {
        type            empty;
    }


Ts:

::

    inlet
    {
        type    zeroGradient;
    }

    outlet
    {
        type    zeroGradient;
    }

    wall_upper
    {
        type	zeroGradient;
    }

    wall_lower
    {
        type	zeroGradient;
    }

    wall_lower_heated
    {
        type groovyBC;
        valueExpression "Tf";
        variables "Tf@wall_upper=Tf;";
        value $internalField;
    }

    frontAndBack
    {
        type            empty;
    }

Tf:

::

    inlet
    {
        type        fixedValue;
        value		uniform 300;
    }

    outlet
    {
        type        zeroGradient;
    }

    wall_upper
    {
        type		zeroGradient;
    }

    wall_lower
    {
        type		zeroGradient;
    }

    wall_lower_heated
    {
        type		fixedGradient;
        gradient	uniform -447.28;
    }

    frontAndBack
    {
        type        empty;
    }

U:

::
    
    inlet
    {
        type            hyperbolicVelocityZhao;
        n               (1 0 0);
        y               (0 1 0);
        maxValue        0.70;
            h			    0.045;
            k       		2.70e-7;
        value           uniform (0 0 0); // Dummy for paraFoam
    }

    outlet
    {
        type        zeroGradient;
    }

    wall_upper
    {
        type        fixedValue;
        value		uniform (0 0 0);
    }

    wall_lower
    {
        type        fixedValue;
        value		uniform (0 0 0);
    }

    wall_lower_heated
    {
        type        fixedValue;
            value		uniform (0 0 0);
    }

    frontAndBack
    {
        type        empty;
    }  
    

porosityProperties:

::

    porosity1
    {
            type            DarcyForchheimerCoefficients;
            active          yes;
            cellZone        porosity;

            DarcyForchheimerCoefficientsCoeffs
            {

                k k [0 2 0 0 0 0 0] 2.70e-7;
                f f [0 0 0 0 0 0 0] 0.097;

                d d [0 -2 0 0 0 0 0] (1.0 1.0 1.0);
                f f [0 -1 0 0 0 0 0] (1.0 1.0 1.0);

                coordinateSystem
                {
                    type    cartesian;
                    origin  (0 0 0);
                    coordinateRotation
                    {
                        type    axesRotation;
                        e1      (1 0 0);
                        e2      (0 1 0);
                    }
                }

            }
    }


transportProperties

::

    transportModel  Newtonian;

    nu             nu [ 0 2 -1 0 0 0 0 ] 1.562e-5;

    // For Darcy-Forchheimer BC:

    k              k  [ 0 2 0 0 0 0 0 ] 2.70e-7;
    
    f              f  [ 0 0 0 0 0 0 0 ] 0.097;

    // For temperature equation:

    eps            eps [ 0 0 0 0 0 0 0 ] 0.9726;

    kfe            kfe [ 1 1 -3 -1 0 0 0 ] 0.0256;

    kse            kse [ 1 1 -3 -1 0 0 0 ] 2.48;

    rhof           rhof [ 1 -3 0 0 0 0 0 ] 1.184;

    cpf            cpf [ 0 2 -2 -1 0 0 0 ] 1007.0;

    hsf	           hsf [ 1 0 -3 -1 0 0 0 ] 1101.412;

    asf		       asf [ 0 -1 0 0 0 0 0 ] 415.4121;

controlDict:

::

    application     porousSimpleTwoTemperaturesZhaoFoam;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         500000;

    deltaT          100;

    writeControl    runTime;

    writeInterval   10000;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;

    libs 
    (
            "ldarcyForchheimerCoefficients.so"
            "libhyperbolicVelocityZhao.so"
            "libmyfiniteVolume.so"
            "libgroovyBC.so"
    );

sampleDict

::

    sets
    (
            lineX1
            {
                type        midPoint;
                axis        distance;

                //- cavity. Slightly perturbed so not to align with face or edge.
                start       (0.0057 0.000 0.000);
                end         (0.1197 0.000 0.000);
            }

    );


Use Momentum Equation of Calmidi and add thermal dispersion
-----------------------------------------------------------

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTwoTemperaturesZhaoFoam simpleTwoTemperaturesCalmidiFoam
    $ cd simpleTwoTemperaturesCalmidiFoam

    $ mv simpleTwoTemperaturesZhaoFoam.C simpleTwoTemperaturesCalmidiFoam.C
    $ sed -i s/simpleTwoTemperaturesZhaoFoam/simpleTwoTemperaturesCalmidiFoam/g Make/files

    $ createFields.H

            Info<< "Reading diffusion conductivity  kd\n" << endl;

            dimensionedScalar kd
            (
                transportProperties.lookup("kd")
            );

    $ nano UEqn.H

            tmp<fvVectorMatrix> UEqn
            (
                (1.0/sqr(eps))*fvm::div(phi, U)
                - fvm::laplacian(nu/eps, U)
                ==
                fvOptions(U)
            );

    $ wclean
    $ wmake

2) Edit the sub-directory

::

    $ mv porousSimpleTwoTemperaturesZhaoFoam porousSimpleTwoTemperaturesCalmidiFoam
    $ cd porousSimpleTwoTemperaturesCalmidiFoam
    $ sed -i s/porousSimpleTwoTemperaturesZhaoFoam/porousSimpleTwoTemperaturesCalmidiFoam/g Make/files
    $ mv porousSimpleTwoTemperaturesZhaoFoam.C porousSimpleTwoTemperaturesCalmidiFoam.C

    $ nano UEqn.H

            tmp<fvVectorMatrix> UEqn
            (
                (1.0/sqr(eps))*fvm::div(phi, U)
                - fvm::laplacian(nu/eps, U)
                ==
                fvOptions(U)
            );

    $ nano TEqn.H

            // Tf treated as implicit in TfEqn, Ts treated as explicit
            fvScalarMatrix TfEqn
            (
                    rhof*cpf*fvm::div(phi,Tf)
                ==
                    fvm::laplacian((kfe+kd),Tf,"laplacian(DT,T)")
                    - fvm::Sp(hsf*asf,Tf)
                    + hsf*asf*Ts
            );
            TfEqn.relax();
            TfEqn.solve();

    $ wclean
    $ wmake


Profile of Vafai
----------------

::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/derived
    $ cp -rf hyperbolicVelocityCalculated hyperbolicVelocityVafai


::

    $ cd hyperbolicVelocityVafai
    $ wclean libso
    $ mv hyperbolicVelocityCalculatedFvPatchVectorField.C hyperbolicVelocityVafaiFvPatchVectorField.C
    $ mv hyperbolicVelocityCalculatedFvPatchVectorField.H hyperbolicVelocityVafaiFvPatchVectorField.H
    $ sed -i s/hyperbolicVelocityCalculated/hyperbolicVelocityVafai/g hyperbolicVelocityVafaiFvPatchVectorField.C
    $ sed -i s/hyperbolicVelocityCalculated/hyperbolicVelocityVafai/g hyperbolicVelocityVafaiFvPatchVectorField.H
    $ sed -i s/hyperbolicVelocityCalculated/hyperbolicVelocityVafai/g Make/files


::

    $ nano hyperbolicVelocityVafaiFvPatchVectorField.H

    //- Peak velocity magnitude
    scalar maxValue_;

    //- Flow direction
    vector n_;

    //- Direction of the y-coordinate
    vector y_;

    //- Forchheimer
            scalar f_;
    
    //- Permeability
            scalar k_;
    
    //- Porosity
    scalar epsilon_;

    //- Plate height
    scalar height_;

    //- Fluid kinematic viscosity
    scalar kinVisc_;


            //- Return max value
    scalar& maxValue()
    {
        return maxValue_;
    }

    //- Return flow direction
    vector& n()
    {
        return n_;
    }

    //- Return y direction
    vector& y()
    {
        return y_;
    }
    
    //- Forchheimer
    scalar& f()
    {
        return f_;
    }

    //- Forchheimer
    scalar& k()
    {
        return k_;
    }
    
    //- Porosity
    scalar& epsilon()
    {
        return epsilon_;
    }
    
    //- height
    scalar& height()
    {
        return height_;
    }

    //- height
    scalar& kinVisc()
    {
        return kinVisc_;
    }

::

	$ nano hyperbolicVelocityVafaiFvPatchVectorField.C

		hyperbolicVelocityVafaiFvPatchVectorField::hyperbolicVelocityVafaiFvPatchVectorField
		(
			const fvPatch& p,
			const DimensionedField<vector, volMesh>& iF
		)
		:
			fixedValueFvPatchVectorField(p, iF),
			maxValue_(0),
			n_(1, 0, 0),
			y_(0, 1, 0),
			f_(1),
			k_(1),
			epsilon_(1),
			height_(1),
			kinVisc_(1)
		{}


		hyperbolicVelocityVafaiFvPatchVectorField::hyperbolicVelocityVafaiFvPatchVectorField
		(
			const hyperbolicVelocityVafaiFvPatchVectorField& ptf,
			const fvPatch& p,
			const DimensionedField<vector, volMesh>& iF,
			const fvPatchFieldMapper& mapper
		)
		:
			fixedValueFvPatchVectorField(ptf, p, iF, mapper),
			maxValue_(ptf.maxValue_),
			n_(ptf.n_),
			y_(ptf.y_),
			f_(ptf.f_),
			k_(ptf.k_),
			epsilon_(ptf.epsilon_),
			height_(ptf.height_),
			kinVisc_(ptf.kinVisc_)
		{}


		hyperbolicVelocityVafaiFvPatchVectorField::hyperbolicVelocityVafaiFvPatchVectorField
		(
			const fvPatch& p,
			const DimensionedField<vector, volMesh>& iF,
			const dictionary& dict
		)
		:
			fixedValueFvPatchVectorField(p, iF),
			maxValue_(readScalar(dict.lookup("maxValue"))),
			n_(dict.lookup("n")),
			y_(dict.lookup("y")),
			f_(readScalar(dict.lookup("f"))),
			k_(readScalar(dict.lookup("k"))),    
			epsilon_(readScalar(dict.lookup("epsilon"))),
			height_(readScalar(dict.lookup("height"))),
			kinVisc_(readScalar(dict.lookup("kinVisc")))
		{
			if (mag(n_) < SMALL || mag(y_) < SMALL)
			{
				FatalErrorIn("hyperbolicVelocityVafaiFvPatchVectorField(dict)")
				    << "n or y given with zero size not correct"
				    << abort(FatalError);
			}

			n_ /= mag(n_);
			y_ /= mag(y_);

			evaluate();
		}


		hyperbolicVelocityVafaiFvPatchVectorField::hyperbolicVelocityVafaiFvPatchVectorField
		(
			const hyperbolicVelocityVafaiFvPatchVectorField& fcvpvf,
			const DimensionedField<vector, volMesh>& iF
		)
		:
			fixedValueFvPatchVectorField(fcvpvf, iF),
			maxValue_(fcvpvf.maxValue_),
			n_(fcvpvf.n_),
			y_(fcvpvf.y_),
			f_(fcvpvf.f_),
			k_(fcvpvf.k_),    
			epsilon_(fcvpvf.epsilon_),
			height_(fcvpvf.height_),
			kinVisc_(fcvpvf.kinVisc_)
		{}


		// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

		void hyperbolicVelocityVafaiFvPatchVectorField::updateCoeffs()
		{
			if (updated())
			{
				return;
			}

			// Get range and orientation
			boundBox bb(patch().patch().localPoints(), true);

			vector ctr = 0.5*(bb.max() + bb.min());

			const vectorField& c = patch().Cf();

			// Calculate local 1-D coordinate for the parabolic profile
			scalarField coord = 2*((c - ctr) & y_)/((bb.max() - bb.min()) & y_);
			   
			const double Lambda = ((pow(epsilon_, (3.0/2.0)))*f_*maxValue_*height_) / kinVisc_;
			const double Darcy = k_ / (sqr(height_)*epsilon_);
			const double ACoeff = (2.0 / 3.0)*Lambda*(pow(Darcy, (-1.0/2.0)));
			const double BCoeff = (pow(Darcy, (-1.0))) + (4.0/3.0)*Lambda*pow(Darcy, (-1.0/2.0));
			const double DCoeff = (sqrt(ACoeff + BCoeff))/2.0; 
			const double Value = sqrt(ACoeff/(ACoeff + BCoeff));
			const double CCoeff = -(1.0 / DCoeff)*acosh(1.0/Value) - 1.0;

			vectorField::operator=(n_*maxValue_*(1.0 - ((ACoeff+BCoeff)/ACoeff)*sqr(  1.0 / (cosh (DCoeff*  (mag(coord) + CCoeff) ) )   )));

			Info << "Permeability" << k_ << endl;
			Info << "Forchheimer" << f_ << endl;
			Info << "Lambda" << Lambda << endl;
			Info << "Darcy" << Darcy << endl;
			Info << "ACoeff" << ACoeff << endl;
			Info << "BCoeff" << BCoeff << endl;
			Info << "DCoeff" << DCoeff << endl;
			Info << "CCoeff" << CCoeff << endl;
		}


		// Write
		void hyperbolicVelocityVafaiFvPatchVectorField::write(Ostream& os) const
		{
			fvPatchVectorField::write(os);
			os.writeKeyword("maxValue")
				<< maxValue_ << token::END_STATEMENT << nl;
			os.writeKeyword("n")
				<< n_ << token::END_STATEMENT << nl;
			os.writeKeyword("y")
				<< y_ << token::END_STATEMENT << nl;
			os.writeKeyword("f")
				<< f_ << token::END_STATEMENT << nl;
			os.writeKeyword("k")
				<< k_ << token::END_STATEMENT << nl;        
			os.writeKeyword("epsilon")
				<< epsilon_ << token::END_STATEMENT << nl;
			os.writeKeyword("height")
				<< height_ << token::END_STATEMENT << nl;
			os.writeKeyword("kinVisc")
				<< kinVisc_ << token::END_STATEMENT << nl;
			writeEntry("value", os);
		}



Recompile:

::

    $ wclean libso
    $ wmake libso


Test Calmidi case
-----------------

Copy case:

::

    $ cp -rf porous_037_porousSimpleTwoTemperaturesZhaoProfile porous_038_porousSimpleTwoTemperaturesCalmidi

U:

::
    
    inlet
    {
        type            hyperbolicVelocityVafai;
        n               (1 0 0);
        y               (0 1 0);
        f        		0.085;
        k       		1.80e-7;
        maxValue        0.61;
        h			    0.045;
        value           uniform (0 0 0); // Dummy for paraFoam
    }

    outlet
    {
        type        zeroGradient;
    }

    wall_upper
    {
        type        fixedValue;
        value	uniform (0 0 0);
    }

    wall_lower
    {
        type     fixedValue;
        value		uniform (0 0 0);
    }

    wall_lower_heated
    {
        type        fixedValue;
		value		uniform (0 0 0);
    }

    frontAndBack
    {
        type        empty;
    }  
    
Tf:

::

    inlet
    {
        type            fixedValue;
        value		uniform 300;
    }

    outlet
    {
        type            zeroGradient;
    }

    wall_upper
    {
        type		zeroGradient;
    }

    wall_lower
    {
        type		zeroGradient;
    }

    wall_lower_heated
    {
        type		fixedGradient;
        gradient	uniform 445.96;
    }

    frontAndBack
    {
        type            empty;
    }

Ts:

::

    inlet
    {
        type    zeroGradient;
    }

    outlet
    {
        type    zeroGradient;
    }
    
    wall_upper
    {
        type	zeroGradient;
    }

    wall_lower
    {
        type	zeroGradient;
    }

    wall_lower_heated
    {
        type groovyBC;
        valueExpression "Tf";
        variables "Tf@wall_upper=Tf;";
        value $internalField;
    }

    frontAndBack
    {
        type            empty;
    }

p:

::

    inlet
    {
        type            darcyForchheimerCoefficients;
        value           uniform 0;
    }

    outlet
    {
        type            fixedMean;
        meanValue       0.0;
        value           uniform 0;
    }

    wall_upper
    {
        type            zeroGradient;
    }

    wall_lower_heated
    {
        type            zeroGradient;
    }

    wall_lower
    {
        type            zeroGradient;
    }

    frontAndBack
    {
        type            empty;
    }

porosityProperties:

::

    DarcyForchheimerCoefficientsCoeffs
    {

        k k [0 2 0 0 0 0 0] 1.80e-7;
        f f [0 0 0 0 0 0 0] 0.085;

        d d [0 -2 0 0 0 0 0] (1.0 1.0 1.0);
        e e [0 -1 0 0 0 0 0] (1.0 1.0 1.0);

        coordinateSystem
        {
            type    cartesian;
            origin  (0 0 0);
            coordinateRotation
            {
                type    axesRotation;
                e1      (1 0 0);
                e2      (0 1 0);
            }
        }

    }

transportProperties

::

    transportModel  Newtonian;

    nu             nu [ 0 2 -1 0 0 0 0 ] 1.608e-5;

    // For Darcy-Forchheimer BC:

    k              k  [ 0 2 0 0 0 0 0 ] 1.80e-7;
    
    f              f  [ 0 0 0 0 0 0 0 ] 0.085;

    // For temperature equation:

    eps            eps [ 0 0 0 0 0 0 0 ] 0.9726;

    kfe            kfe [ 1 1 -3 -1 0 0 0 ] 0.0237;

    kse            kse [ 1 1 -3 -1 0 0 0 ] 6.46;

    rhof           rhof [ 1 -3 0 0 0 0 0 ] 1.164;

    cpf            cpf [ 0 2 -2 -1 0 0 0 ] 1007.0;

    hsf	           hsf [ 1 0 -3 -1 0 0 0 ] 104.08;

    asf		   	   asf [ 0 -1 0 0 0 0 0 ] 611.25;

controlDict

::

    application     porousSimpleTwoTemperaturesCalmidiFoam;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         500000;

    deltaT          100;

    writeControl    runTime;

    writeInterval   10000;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;

    libs 
    (
            "ldarcyForchheimerCoefficients.so"
            "libhyperbolicVelocityVafai.so"
            "libmyfiniteVolume.so"
            "libgroovyBC.so"
    );


Copy in residuals:

::

    residuals
    {
            fields ( p U Tf Ts );

            #includeEtc "caseDicts/postProcessing/numerical/residuals.cfg"
    }

sampleDict:

::

    lineX1
    {
        type        midPoint;
        axis        distance;

        //- cavity. Slightly perturbed so not to align with face or edge - wall.
        start       (0.000 0.000 0.000);
        end         (0.114 0.000 0.000);
    }
    lineY1
    {
        type        midPoint;
        axis        distance;

        //- cavity. Slightly perturbed so not to align with face or edge - inlet.
        start       (0.000 0.000 0.000);
        end         (0.000 0.045 0.000);
    }
    lineY2
    {
        type        midPoint;
        axis        distance;

        //- cavity. Slightly perturbed so not to align with face or edge - centre.
        start       (0.057 0.000 0.000);
        end         (0.057 0.045 0.000);
    }
    lineY3
    {
        type        midPoint;
        axis        distance;

        //- cavity. Slightly perturbed so not to align with face or edge - outlet.
        start       (0.114 0.000 0.000);
        end         (0.114 0.045 0.000);
    }

controlDict

::

    application     porousSimpleTwoTemperaturesCalmidiFoam;

    startFrom       startTime;

    startTime       0;

    stopAt          endTime;

    endTime         500000;

    deltaT          100;

    writeControl    runTime;

    writeInterval   10000;

    purgeWrite      0;

    writeFormat     ascii;

    writePrecision  6;

    writeCompression off;

    timeFormat      general;

    timePrecision   6;

    runTimeModifiable true;

    libs 
    (
            "ldarcyForchheimerCoefficients.so"
            "libhyperbolicVelocityVafai.so"
            "libmyfiniteVolume.so"
            "libgroovyBC.so"
    );

    functions
    {
            residuals
            {
                type            residuals;
                functionObjectLibs ("libutilityFunctionObjects.so");
                enabled         true;
                outputControl   timeStep;
                outputInterval  1;

                fields
                (
                    Ux
                    Uy
                    Uz
                    p
                    Tf
                        Ts
                );
            }

            probes
            {
                type		      probes;
                functionObjectLibs    ("libsampling.so");
                enabled		      true;
                outputControl	      timeStep;
                outputInterval	      1;

                fields
                (
                    Ux
                    Uy
                    Uz
                    p
                    Tf
                        Ts
                );
                    
                probeLocations
                (
                    (0.0114 0.0225 0.0) // Probe 1 centreline
                    (0.0570 0.0225 0.0) // Probe 2 centreline
                    (0.1026 0.0225 0.0) // Probe 3 centreline
                    (0.0114 0.0 0.0) // Probe 4 wall 
                    (0.0570 0.0 0.0) // Probe 5 wall
                    (0.1026 0.0 0.0) // Probe 6 wall
                );
            }
    }

Monitor residuals:

::

    $ foamMonitor -l postProcessing/residuals/0/residuals.dattree
    $ foamMonitor -l postProcessing/probes/0/p
    $ foamMonitor -l postProcessing/probes/0/Tf
    $ foamMonitor -l postProcessing/probes/0/Ts

    $ porousSimpleTwoTemperaturesCalmidiFoam
    $ rm -rf 0.* [1-9]* postProcessing processor* *.foam logs *.png

    $ foamCalc mag U

    $ gnuplot -persist plot_velocity_probes


Introduce transient
-------------------

::

    $ cp -rf porous_038_porousSimpleTwoTemperaturesCalmidi porous_039_porousSimpleTwoTemperaturesCalmidiTransient

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTwoTemperaturesCalmidiFoam simpleTwoTemperaturesCalmidiTransientFoam
    $ cd simpleTwoTemperaturesCalmidiTransientFoam

    $ mv simpleTwoTemperaturesCalmidiFoam.C simpleTwoTemperaturesCalmidiTransientFoam.C
    $ sed -i s/simpleTwoTemperaturesCalmidiFoam/simpleTwoTemperaturesCalmidiTransientFoam/g Make/files


::

    $ mv porousSimpleTwoTemperaturesCalmidiFoam porousSimpleTwoTemperaturesCalmidiTransientFoam
    $ cd porousSimpleTwoTemperaturesCalmidiTransientFoam
    $ sed -i s/porousSimpleTwoTemperaturesCalmidiFoam/porousSimpleTwoTemperaturesCalmidiTransientFoam/g Make/files
    $ mv porousSimpleTwoTemperaturesCalmidiFoam.C porousSimpleTwoTemperaturesCalmidiTransientFoam.C

    $ kate TEqn.H

    fvScalarMatrix TfEqn
    (
          eps*(rhos*cpf)*fvm::ddt(Tf)
        + rhof*cpf*fvm::div(phi,Tf)
        ==
          fvm::laplacian((kfe+kd),Tf,"laplacian(DT,T)")
        - fvm::Sp(hsf*asf,Tf)
        + hsf*asf*Ts
    );
    
    TfEqn.relax();
    TfEqn.solve();

    // Ts treated an implicit in TsEqn, Tf treated as explicit
    fvScalarMatrix TsEqn
    (
        (1.0 - eps)*(rhos*cps)*fvm::ddt(Ts)
        ==
          fvm::laplacian(kse,Ts,"laplacian(DT,T)")
        - fvm::Sp(hsf*asf,Ts)
        + hsf*asf*Tf
    );
    TsEqn.relax();
    TsEqn.solve();

    $ wclean
    $ wmake

Add to createFields.H

::

    Info<< "Reading fluid density rhos\n" << endl;

    dimensionedScalar rhos
    (
        transportProperties.lookup("rhos")
    );

    Info<< "Reading fluid specific heat capacity cps\n" << endl;

    dimensionedScalar cps
    (
        transportProperties.lookup("cps")
    );

    $ wclean
    $ wmake

script:

::

    gnuplot plot_p_probes
    gnuplot plot_residuals
    gnuplot plot_Tf_probes
    gnuplot plot_Ts_probes
    gnuplot plot_velocity_probes

    $ chmod 777

    $ rm -rf 0.* [1-9]* postProcessing processor* *.foam logs *.png


    $ porousSimpleTwoTemperaturesCalmidiTransientFoam > log &


    $ sampleDict




Influence of the grid
---------------------

Velocity profile vs Vafai

::

    $ cp -rf porous_039_porousSimpleTwoTemperaturesCalmidiTransient porous_040_porousSimpleTwoTemperaturesCalmidiTransientMesh
    $ rm -rf 0.* [1-9]* postProcessing processor* *.foam logs *.png

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleTwoTemperaturesCalmidiFoam simpleTwoTemperaturesCalmidiTransientFoam
    $ cd simpleTwoTemperaturesCalmidiTransientFoam

    $ mv simpleTwoTemperaturesCalmidiFoam.C simpleTwoTemperaturesCalmidiTransientFoam.C
    $ sed -i s/simpleTwoTemperaturesCalmidiFoam/simpleTwoTemperaturesCalmidiTransientFoam/g Make/files




Message: Number of cells/points in mesh and field don't match - something in 0 folder e.g. magU doesn;t match the grid.


Changed inlet pressure condition to fixedFluxPressure
-----------------------------------------------------

So that pressure gradient is computed based on velocity - didn't work because updateCoeffs for snGradp must be done before updateCoeffs or evaluate.

Tried porousSimpleFoam
----------------------

Almost matches maximum. Slighly less diffuse than Vafai.


Now modify porousSimpleFoam so that diffusion is included in laplacian
----------------------------------------------------------------------

Use Momentum Equation of Vafai
------------------------------

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleFoam simpleEpsilonFoam
    $ cd simpleEpsilonFoam

    $ mv simpleFoam.C simpleEpsilonFoam.C
    $ sed -i s/simpleFoam/simpleEpsilonFoam/g Make/files

    $ nano createFields.H

            Info<< "Reading transportProperties\n" << endl;

            IOdictionary transportProperties
            (
                IOobject
            (
                    "transportProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
            )
            );

            Info<< "Reading kinematic viscosity nu\n" << endl;

            dimensionedScalar nu
            (
                transportProperties.lookup("nu")
            );

            Info<< "Reading fluid viscosity eps\n" << endl;

            dimensionedScalar eps
            (
                transportProperties.lookup("eps")
            );

::

    $ mv porousSimpleFoam porousSimpleEpsilonFoam
    $ cd porousSimpleEpsilonFoam
    $ sed -i s/porousSimpleFoam/porousSimpleEpsilonFoam/g Make/files
    $ mv porousSimpleFoam.C porousSimpleEpsilonFoam.C

    $ nano UEqn.H

            tmp<fvVectorMatrix> UEqn
            (
                fvm::div(phi, U)
                - fvm::laplacian(nu/eps, U)
                ==
                fvOptions(U)
            );

    $ wclean
    $ wmake

Use Momentum Equation of Vafai - convection term over epsilon squared
---------------------------------------------------------------------

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -rf simpleEpsilonFoam simpleEpsilonSquaredFoam
    $ cd simpleEpsilonSquaredFoam

    $ mv simpleEpsilonFoam.C simpleEpsilonSquaredFoam.C
    $ sed -i s/simpleEpsilonFoam/simpleEpsilonSquaredFoam/g Make/files

    $ wclean
    $ wmake


::

    $ mv porousSimpleEpsilonFoam porousSimpleEpsilonSquaredFoam
    $ cd porousSimpleEpsilonSquaredFoam
    $ sed -i s/porousSimpleEpsilonFoam/porousSimpleEpsilonSquaredFoam/g Make/files
    $ mv porousSimpleEpsilonFoam.C porousSimpleEpsilonSquaredFoam.C

    $ nano UEqn.H

            tmp<fvVectorMatrix> UEqn
            (
                fvm::div(phi/sqr(eps), U)
                - fvm::laplacian(nu/eps, U)
                ==
                fvOptions(U)
            );

    $ wclean
    $ wmake

Vafai includes epsilon next to Forchheimer term
-----------------------------------------------

::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/cfdTools/general/porosityModel
    $ cp -rf DarcyForchheimerCoefficients DarcyForchheimerCoefficientsEpsilon
    $ cd DarcyForchheimerCoefficientsEpsilon

    $ mv DarcyForchheimerCoefficients.C DarcyForchheimerCoefficientsEpsilon.C
    $ mv DarcyForchheimerCoefficients.H DarcyForchheimerCoefficientsEpsilon.H
    $ mv DarcyForchheimerCoefficientsTemplates.C DarcyForchheimerCoefficientsEpsilonTemplates.C
    
::

    $ sed -i s/DarcyForchheimerCoefficients/DarcyForchheimerCoefficientsEpsilon/g DarcyForchheimerCoefficientsEpsilon.C
    $ sed -i s/DarcyForchheimerCoefficients/DarcyForchheimerCoefficientsEpsilon/g DarcyForchheimerCoefficientsEpsilon.H
    $ sed -i s/DarcyForchheimerCoefficients/DarcyForchheimerCoefficientsEpsilon/g DarcyForchheimerCoefficientsEpsilonTemplates.C










::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/Make
    $ gedit files

    general = cfdTools/general
    porosity = $(general)/porosityModel
    $(porosity)/DarcyForchheimerCoefficientsEpsilon/DarcyForchheimerCoefficientsEpsilon.C
    LIB = $(FOAM_USER_LIBBIN)/libmyfiniteVolume

    $ gedit options

            EXE_INC = \
    -I$(LIB_SRC)/triSurface/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/finiteVolume/lnInclude
    LIB_LIBS = \
    -lOpenFOAM \
    -ltriSurface \
    -lmeshTools


::

    $ cd ..
    $ wclean libso
    $ wmake libso

::

    porosity1
    {
        type            DarcyForchheimerCoefficientsEpsilon;
        active          yes;
        cellZone        porosity;

        DarcyForchheimerCoefficientsEpsilonCoeffs
        {

            k k [0 2 0 0 0 0 0] 1.80e-7;
            f f [0 0 0 0 0 0 0] 0.077503;

            d d [0 -2 0 0 0 0 0] (1.0 1.0 1.0);
            e e [0 -1 0 0 0 0 0] (1.0 1.0 1.0);

            coordinateSystem
            {
                type    cartesian;
                origin  (0 0 0);
                coordinateRotation
                {
                    type    axesRotation;
                    e1      (1 0 0);
                    e2      (0 1 0);
                }
            }

        }
    }





::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/
    $ wclean libso



11/09/2017: Add temperature to solver (local equilibirum assumption):
---------------------------------------------------------------------

simpleEpsilonSquaredTemperatureFoam
porousSimpleEpsilonSquaredTemperatureFoam

::

    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $WM_PROJECT_USER_DIR/applications/solvers/simpleEpsilonSquaredFoam simpleEpsilonSquaredTemperatureFoam
    $ cd simpleEpsilonSquaredTemperatureFoam
    $ mv simpleEpsilonSquaredFoam.C simpleEpsilonSquaredTemperatureFoam.C
    $ sed -i s/simpleEpsilonSquaredFoam/simpleEpsilonSquaredTemperatureFoam/g Make/files

::

    $ mv porousSimpleEpsilonSquaredFoam porousSimpleEpsilonSquaredTemperatureFoam
    $ cd porousSimpleEpsilonSquaredTemperatureFoam
    $ sed -i s/porousSimpleEpsilonSquaredFoam/porousSimpleEpsilonSquaredTemperatureFoam/g Make/files
    $ mv porousSimpleEpsilonSquaredFoam.C porousSimpleEpsilonSquaredTemperatureFoam.C

::

    $ rm -rf 0.* [1-9]* postProcessing processor* *.foam logs *.png

Add temperature field to createFields.H

::

    $ kate createFields.H

        Info<< "Reading field T\n" << endl;
        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        #include "createPhi.H"

        Info<< "Reading effective fluid conductivity kfe\n" << endl;

        dimensionedScalar kfe
        (
            transportProperties.lookup("kfe")
        );
        
        Info<< "Reading effective solid conductivity kse\n" << endl;

        dimensionedScalar kse
        (
            transportProperties.lookup("kse")
        );    
        
        Info<< "Reading fluid density rhof\n" << endl;

        dimensionedScalar rhof
        (
            transportProperties.lookup("rhof")
        );
        
        Info<< "Reading fluid specific heat capacity cpf\n" << endl;

        dimensionedScalar cpf
        (
            transportProperties.lookup("cpf")
        );

Add TEqn.H    
    
::

    $ kate TEqn.H

        fvScalarMatrix TEqn
        (
            fvm::div(phi,T)
            ==
            fvm::laplacian((kfe+kse)/(rhof*cpf*eps),T,"laplacian(DT,T)")
        );
        TEqn.relax();
        TEqn.solve();
    
sampleDict modification
    
::
    
    $ kate sampleDict

        fields
        (
            p
            U
            T
        );

controlDict modification        

::

    $ kate controlDict

        residuals
        {
            type            residuals;
            functionObjectLibs ("libutilityFunctionObjects.so");
            enabled         true;
            outputControl   timeStep;
            outputInterval  1;

            fields
            (
                U
                p
                T
            );
        }
        
        probes
        {
    
            type		      probes;
            functionObjectLibs    ("libsampling.so");
            enabled		      true;
            outputControl	      timeStep;
            outputInterval	      1;
            fields
            (
                U
                p
                T
            );
            
            probeLocations
            (
                (0.0049 0.0225 0.0) // Probe 1 centreline
                (0.0189 0.0225 0.0) // Probe 2 centreline
                (0.0386 0.0225 0.0) // Probe 3 centreline
                (0.0816 0.0225 0.0) // Probe 4 centreline
                (0.1054 0.0225 0.0) // Probe 5 centreline
                (0.0049 0.0 0.0) // Probe 1 wall 
                (0.0189 0.0 0.0) // Probe 2 wall
                (0.0386 0.0 0.0) // Probe 3 wall
                (0.0816 0.0 0.0) // Probe 4 wall
                (0.1054 0.0 0.0) // Probe 5 wall
            );
        }         
        
Temperature doesn't change

controlDict

::
        
    ddtSchemes
    {
        default         steadyState;
    }
    
    deltaT          1;
    
fvSchemes    

::
      
    div(phi,T) bounded Gauss upwind;      
      
porousSimpleEpsilonSquaredTemperatureFoam.C - include TEqn.H   

::

    while (simple.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            // Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }

            turbulence->correct();
            
            #include "TEqn.H"
            
            runTime.write();

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
        
fvSolution

::

    T
    {
        solver          PBiCG;
        preconditioner  DILU;
        tolerance       1e-06;
        relTol          0;
    }          
        


28/10/2017: Change viscosity name:
----------------------------------

simpleEpsilonSquaredTemperatureConductivityFoam
porousSimpleEpsilonSquaredTemperatureConductivityFoam

::

    $ two
    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $WM_PROJECT_USER_DIR/applications/solvers/simpleEpsilonSquaredTemperatureFoam simpleEpsilonSquaredTemperatureConductivityFoam
    $ cd simpleEpsilonSquaredTemperatureConductivityFoam
    $ mv simpleEpsilonSquaredTemperatureFoam.C simpleEpsilonSquaredTemperatureConductivityFoam.C
    $ sed -i s/simpleEpsilonSquaredTemperatureFoam/simpleEpsilonSquaredTemperatureConductivityFoam/g Make/files

    $ wclean
    $ wmake

::

    $ mv porousSimpleEpsilonSquaredTemperatureFoam porousSimpleEpsilonSquaredTemperatureConductivityFoam
    $ cd porousSimpleEpsilonSquaredTemperatureConductivityFoam
    $ sed -i s/porousSimpleEpsilonSquaredTemperatureFoam/porousSimpleEpsilonSquaredTemperatureConductivityFoam/g Make/files
    $ mv porousSimpleEpsilonSquaredTemperatureFoam.C porousSimpleEpsilonSquaredTemperatureConductivityFoam.C

    $ wclean
    $ wmake

Find 

::

    $ find -name '*.H*' | xargs grep -n '.*nu.*'
    
    
::

    $ cd $WM_PROJECT_USER_DIR/src/finiteVolume/cfdTools/general/porosityModel
    $ cp -rf DarcyBrinkmannForchheimer DarcyBrinkmannForchheimerConductivity
    $ cd DarcyBrinkmannForchheimerConductivity

    $ mv DarcyForchheimer.C DarcyBrinkmannForchheimerConductivity.C
    $ mv DarcyForchheimer.H DarcyBrinkmannForchheimerConductivity.H
    $ mv DarcyForchheimerTemplates.C DarcyBrinkmannForchheimerConductivityTemplates.C
    
::

    $ sed -i s/DarcyForchheimer/DarcyBrinkmannForchheimerConductivity/g DarcyBrinkmannForchheimerConductivity.C
    $ sed -i s/DarcyForchheimer/DarcyBrinkmannForchheimerConductivity/g DarcyBrinkmannForchheimerConductivity.H
    $ sed -i s/DarcyForchheimer/DarcyBrinkmannForchheimerConductivity/g DarcyBrinkmannForchheimerConductivityTemplates.C    
    
    changed nu to nulaminar so it is read by dictionary
    
::
    
    $ cd Make
    general = cfdTools/general
    porosity = $(general)/porosityModel
    $(porosity)/DarcyBrinkmannForchheimerConductivity/DarcyBrinkmannForchheimerConductivity.C
    LIB = $(FOAM_USER_LIBBIN)/libmyfiniteVolume
    
    $ cd ..
    $ wclean libso
    $ wmake libso    
    
    
    porosity1
    {
        type            DarcyBrinkmannForchheimerConductivity;
        active          yes;
        cellZone        porosity;

        DarcyForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (4.20323662e9 4.20323662e9 4.20323662e9);
            f   f [0 -1 0 0 0 0 0] (30055.40166205 30055.40166205 30055.40166205);

            coordinateSystem
            {
                type    cartesian;
                origin  (0 0 0);
                coordinateRotation
                {
                    type    axesRotation;
                    e1      (1 0 0);
                    e2      (0 1 0);
                }
            }
        }
    }
    
    Added to controlDict:
    
        "libmyfiniteVolume.so"
        
    Added to porosityProperties:
        
        DarcyBrinkmannForchheimerConductivity


31/10/2017: Add dictionary entires dt and w:
--------------------------------------------

simpleEpsilonSquaredTemperatureConductivityWallFoam
porousSimpleEpsilonSquaredTemperatureConductivityWallFoam

::

    $ two
    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $WM_PROJECT_USER_DIR/applications/solvers/simpleEpsilonSquaredTemperatureConductivityFoam simpleEpsilonSquaredTemperatureConductivityWallFoam
    $ cd simpleEpsilonSquaredTemperatureConductivityWallFoam
    $ mv simpleEpsilonSquaredTemperatureConductivityFoam.C simpleEpsilonSquaredTemperatureConductivityWallFoam.C
    $ sed -i s/simpleEpsilonSquaredTemperatureConductivityFoam/simpleEpsilonSquaredTemperatureConductivityWallFoam/g Make/files

    $ wclean
    $ wmake

::

    $ mv porousSimpleEpsilonSquaredTemperatureConductivityFoam porousSimpleEpsilonSquaredTemperatureConductivityWallFoam
    $ cd porousSimpleEpsilonSquaredTemperatureConductivityWallFoam
    $ sed -i s/porousSimpleEpsilonSquaredTemperatureConductivityFoam/porousSimpleEpsilonSquaredTemperatureConductivityWallFoam/g Make/files
    $ mv porousSimpleEpsilonSquaredTemperatureConductivityFoam.C porousSimpleEpsilonSquaredTemperatureConductivityWallFoam.C

    $ wclean
    $ wmake

Find 

    $ find -name '*.H*' | xargs grep -n '.*nu.*'

Add dt and w to createFields.H

::

    $ kate createFields.H

        Info<< "Reading dt\n" << endl;

        dimensionedScalar dt
        (
            transportProperties.lookup("dt")
        );

        Info<< "Reading w\n" << endl;

        dimensionedScalar w
        (
            transportProperties.lookup("w")
        );

Add dt and w to equations in TEqn.H:

::

    kd=d_t*rhof*cpf*U.component(0)*dp*(1.0-exp(-mesh.C().component(1)/(w*dp)));
    
    kd_2=d_t*rhof*cpf*U.component(0)*dp*(1.0-exp(-(H-mesh.C().component(1))/(w*dp)));

Two temperatures
----------------

simpleTwoTemperaturesFoam
porousSimpleTwoTemperaturesFoam

::

    $ two
    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $WM_PROJECT_USER_DIR/applications/solvers/simpleEpsilonSquaredTemperatureConductivityWallFoam simpleTwoTemperaturesFoam
    $ cd simpleTwoTemperaturesFoam
    $ mv simpleEpsilonSquaredTemperatureConductivityWallFoam.C simpleTwoTemperaturesFoam.C
    $ sed -i s/simpleEpsilonSquaredTemperatureConductivityWallFoam/simpleTwoTemperaturesFoam/g Make/files

    $ wclean
    $ wmake

::

    $ mv porousSimpleEpsilonSquaredTemperatureConductivityWallFoam porousSimpleTwoTemperaturesFoam
    $ cd porousSimpleTwoTemperaturesFoam
    $ sed -i s/porousSimpleEpsilonSquaredTemperatureConductivityWallFoam/porousSimpleTwoTemperaturesFoam/g Make/files
    $ mv porousSimpleEpsilonSquaredTemperatureConductivityWallFoam.C porousSimpleTwoTemperaturesFoam.C

    $ wclean
    $ wmake

TEqn.H

::

    kd = (kse+kfe)*0.1*U.component(0)/u_mean;
    
    // update coupled boundaries
    kd.correctBoundaryConditions();

    //record updated kd
    if(runTime.outputTime())
    {
      kd.write();
    }

    // Tf treated as implicit in TfEqn, Ts treated as explicit
    fvScalarMatrix TfEqn
    (
          rhof*cpf*fvm::div(phi,Tf)
        ==
          fvm::laplacian((kfe+kd),Tf,"laplacian(DT,T)")
        - fvm::Sp(hsf*asf,Tf)
        + hsf*asf*Ts
    );

    TfEqn.relax();
    TfEqn.solve();

    // Ts treated an implicit in TsEqn, Tf treated as explicit
    fvScalarMatrix TsEqn
    (
          fvm::laplacian(kse,Ts,"laplacian(DT,T)")
        - fvm::Sp(hsf*asf,Ts)
        + hsf*asf*Tf
    );

    TsEqn.relax();
    TsEqn.solve();

createFields.H

::

    Info<< "Reading interfacial heat transfer coefficient hsf\n" << endl;

    dimensionedScalar hsf
    (
        transportProperties.lookup("hsf")
    );

    Info<< "Reading surface area density asf\n" << endl;

    dimensionedScalar asf
    (
        transportProperties.lookup("asf")
    );


createFields.H

::

    Info<< "Reading field Ts\n" << endl;
    volScalarField Ts
    (
        IOobject
        (
            "Ts",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field Tf\n" << endl;
    volScalarField Tf
    (
        IOobject
        (
            "Tf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

Two temperatures plus dispersion
--------------------------------

simpleTwoTemperaturesDispersionFoam
porousSimpleTwoTemperaturesDispersionFoam


::

    $ two
    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $WM_PROJECT_USER_DIR/applications/solvers/simpleTwoTemperaturesFoam simpleTwoTemperaturesDispersionFoam
    $ cd simpleTwoTemperaturesDispersionFoam 
    $ mv simpleTwoTemperaturesFoam.C simpleTwoTemperaturesDispersionFoam.C
    $ sed -i s/simpleTwoTemperaturesFoam/simpleTwoTemperaturesDispersionFoam/g Make/files

    $ wclean
    $ wmake

::

    $ mv porousSimpleTwoTemperaturesFoam porousSimpleTwoTemperaturesDispersionFoam
    $ cd porousSimpleTwoTemperaturesDispersionFoam
    $ sed -i s/porousSimpleTwoTemperaturesFoam/porousSimpleTwoTemperaturesDispersionFoam/g Make/files
    $ mv porousSimpleTwoTemperaturesFoam.C porousSimpleTwoTemperaturesDispersionFoam.C

    $ wclean
    $ wmake


Two temperatures axi symmetric
------------------------------

simpleTwoTemperaturesDispersionAxiSymmetricFoam
porousSimpleTwoTemperaturesDispersionAxiSymmetricFoam


::

    $ two
    $ cd $WM_PROJECT_USER_DIR/applications/solvers/
    $ cp -r $WM_PROJECT_USER_DIR/applications/solvers/simpleTwoTemperaturesDispersionFoam simpleTwoTemperaturesDispersionAxiSymmetricFoam
    $ cd simpleTwoTemperaturesDispersionAxiSymmetricFoam 
    $ mv simpleTwoTemperaturesDispersionFoam.C simpleTwoTemperaturesDispersionAxiSymmetricFoam.C
    $ sed -i s/simpleTwoTemperaturesDispersionFoam/simpleTwoTemperaturesDispersionAxiSymmetricFoam/g Make/files

    $ wclean
    $ wmake

::

    $ mv porousSimpleTwoTemperaturesDispersionFoam porousSimpleTwoTemperaturesDispersionAxiSymmetricFoam
    $ cd porousSimpleTwoTemperaturesDispersionAxiSymmetricFoam
    $ mv porousSimpleTwoTemperaturesDispersionFoam.C porousSimpleTwoTemperaturesDispersionAxiSymmetricFoam.C
    $ sed -i s/porousSimpleTwoTemperaturesDispersionFoam/porousSimpleTwoTemperaturesDispersionAxiSymmetricFoam/g Make/files
    

    $ wclean
    $ wmake

