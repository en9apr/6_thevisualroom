============
Powder Notes
============

.. contents::
   :local:

Tutorial on twoPhaseEulerFoam
=============================

Drag models: Schiller-Naumann, Wen-Yu and Syamlal-O'Brien

Change to run directory:

::

	$ cd $FOAM_RUN

Copy tutorial files from tutorial directory to new directory:

::

	$ cp -r $FOAM_TUTORIALS/multiphase/twoPhaseEulerFoam/RAS/fluidisedBed $FOAM_RUN


Rename directory:

::

	$ mv fluidisedBed eulerian_001_fluidised_bed

Change directory:

::

	$ cd eulerian_001_fluidised_bed

Mesh
----

Cell size is 30 times particle diameter of solid phase. Cell size cannot be smaller than particle diameter - if it were, this would mean the particle fills the volume and no particle-particle interactions would be possible. This would also violate the maximum solid packing fraction, intended to be 62%. If this limit is exceeded, the probability of collisions between particles would go to infinity.  


::

	$ nano constant/polyMesh/blockMeshDict

		vertices
		(
			(0 0 0)
			(0.28 0 0)
			(0.28 1 0)
			(0 1 0)
			(0 0 0.025)
			(0.28 0 0.025)
			(0.28 1 0.025)
			(0 1 0.025)
		);

		blocks
		(
			hex (0 1 2 3 4 5 6 7) (28 100 1) simpleGrading (1 1 1)
		);


Run blockMesh:

::

	$ blockMesh

Boundary and initial conditions
-------------------------------

Initial conditions
~~~~~~~~~~~~~~~~~~

Create initial volume fraction of particles up to 0.4m

::

	$ system/setFieldsDict

		defaultFieldValues
		(
			volScalarFieldValue alpha.air 1
			volScalarFieldValue alpha.particles 0
		);

		regions
		(
			boxToCell
			{
				box ( 0 0 0 ) ( 0.28 0.4 0.025 );
				fieldValues
				(
				    volScalarFieldValue alpha.air 0.60
				    volScalarFieldValue alpha.particles 0.40
				);
			}
		);

Copy original files (with internalField = 0), as if the (internalField != 0) has not been set.

::

	$ cp 0/alpha.air.org 0/alpha.air
	$ cp 0/alpha.particles.org 0/alpha.particles

Run setFields in order to set the internalField:

::

	$ setFields


Boundary conditions
~~~~~~~~~~~~~~~~~~~

Set the superficial gas velocity to 0.46m/s in 0/U.air:

::

	$ nano 0/U.air

		internalField   uniform (0 0.46 0);

		boundaryField
		{
			inlet
			{
				type               interstitialInletVelocity;
				inletVelocity      uniform (0 0.46 0);
				alpha              alpha.air;
				value              $internalField;
			}

			...

		}


Physical Properties
-------------------

Particle diameter
~~~~~~~~~~~~~~~~~

Change the particle diameter and the maximum solid packing fraction:

::

	$ nano constant/phaseProperties

		phases (particles air);

		particles
		{
			diameterModel constant;
			constantCoeffs
			{
				d               275e-6;
			}

			alphaMax 0.62;
		}


Drag Model
~~~~~~~~~~

Specify drag model (SyamlalOBrien):

::

	$ nano constant/phaseProperties

		drag
		(
			(particles in air)
			{
				type            SyamlalOBrien;
				residualAlpha   1e-6;
				residualRe      1e-3;
				swarmCorrection
				{
				    type        none;
				}
			}
		);

Thermophysical properties
~~~~~~~~~~~~~~~~~~~~~~~~~

Set particle density

::

	$ nano constant/thermophysicalProperties.particles

		equationOfState
			{
				rho         2500;
			}

Turbulence properties
~~~~~~~~~~~~~~~~~~~~~

Set turbulence properties

::

	$ nano constant/turbulenceProperties.particles
	$ nano constant/turbulenceProperties.air


Gravity properties
~~~~~~~~~~~~~~~~~~

::

	$ nano constant/g

Set the endTime to 20s
~~~~~~~~~~~~~~~~~~~~~~

::

	$ nano system/controlDict

		endTime         20;

Parallel processing
~~~~~~~~~~~~~~~~~~~

Create dictionary for parallel processing (scotch - no user input, will minimise the number of processor boundaries) 

::

	numberOfSubdomains 8;

	method          scotch;

	simpleCoeffs
	{
		n               ( 2 1 1 );
		delta           0.0001;
	}

	hierarchicalCoeffs
	{
		n               ( 2 1 1 );
		delta           0.001;
		order           xyz;
	}

	manualCoeffs
	{
		dataFile        "cellDecomposition";
	}

Decompose the case:

::

	$ decomposePar

Run the case:

::

	$ mpirun -np 8 twoPhaseEulerFoam -parallel > log &


Run the case
------------

::

	$ twoPhaseEulerFoam >& log &

Monitor the run

::

	$ tail -f log

::

	$ rm -rf processor* 0.* [1-9]* postProcessing *.foam


Reconstruct the case
--------------------

::

	$ reconstructPar



Created another mesh in Pointwise
---------------------------------
::

    $ $HOME/Pointwise_User/eulerian_002_bin_discharge/pointwise -b script10.glf

::

    $ mkdir eulerian_002_bin_discharge
    
Copied files from $HOME/Pointwise_User/eulerian_002_bin_discharge/mesh/ to current directory

::

    $ cd $FOAM_RUN/eulerian_002_bin_discharge/eulerian_002_bin_discharge


Ran check mesh

::

    $ checkMesh


In the dictionaries:

    epsilon.air, k.air, p, T.air, Theta.particles, U.air :

    inlet=outlet
    outlet=inlet


    All files:

    frontAndBackPlanes=symmetry
    walls=wall


    sed -i s/frontAndBackPlanes/symmetry/g alpha.air.org
    sed -i s/frontAndBackPlanes/symmetry/g alpha.particles.org
    sed -i s/frontAndBackPlanes/symmetry/g epsilon.air
    sed -i s/frontAndBackPlanes/symmetry/g k.air
    sed -i s/frontAndBackPlanes/symmetry/g nut.air
    sed -i s/frontAndBackPlanes/symmetry/g nut.particles
    sed -i s/frontAndBackPlanes/symmetry/g p
    sed -i s/frontAndBackPlanes/symmetry/g T.air
    sed -i s/frontAndBackPlanes/symmetry/g T.particles
    sed -i s/frontAndBackPlanes/symmetry/g Theta.particles
    sed -i s/frontAndBackPlanes/symmetry/g U.air
    sed -i s/frontAndBackPlanes/symmetry/g U.particles

    sed -i s/walls/wall/g alpha.air.org
    sed -i s/walls/wall/g alpha.particles.org
    sed -i s/walls/wall/g epsilon.air
    sed -i s/walls/wall/g k.air
    sed -i s/walls/wall/g nut.air
    sed -i s/walls/wall/g nut.particles
    sed -i s/walls/wall/g p
    sed -i s/walls/wall/g T.air
    sed -i s/walls/wall/g T.particles
    sed -i s/walls/wall/g Theta.particles
    sed -i s/walls/wall/g U.air
    sed -i s/walls/wall/g U.particles


Copy original files (with internalField = 0), as if the (internalField != 0) has not been set.

::

	$ cp 0/alpha.air.org 0/alpha.air
	$ cp 0/alpha.particles.org 0/alpha.particles

Run setFields in order to set the internalField:

::

	$ setFields

Set the endTime to 2s
~~~~~~~~~~~~~~~~~~~~~~

::

	$ nano system/controlDict

		endTime         2;


Run the case
------------

::

	$ twoPhaseEulerFoam > log &




#0  Foam::error::printStack(Foam::Ostream&) at ??:?
#1  Foam::sigFpe::sigHandler(int) at ??:?
#2  ? in "/lib/x86_64-linux-gnu/libc.so.6"
#3  Foam::symGaussSeidelSmoother::smooth(Foam::word const&, Foam::Field<double>&, Foam::lduMatrix const&, Foam::Field<double> const&, Foam::FieldField<Foam::Field, double> const&, Foam::UPtrList<Foam::lduInterfaceField const> const&, unsigned char, int) at ??:?
#4  Foam::symGaussSeidelSmoother::smooth(Foam::Field<double>&, Foam::Field<double> const&, unsigned char, int) const at ??:?
#5  Foam::smoothSolver::solve(Foam::Field<double>&, Foam::Field<double> const&, unsigned char) const at ??:?
#6  Foam::fvMatrix<double>::solveSegregated(Foam::dictionary const&) at ??:?
#7  Foam::fvMatrix<double>::solve(Foam::dictionary const&) at ??:?
#8  Foam::fvMatrix<double>::solve() at ??:?
#9  ? at ??:?
#10  __libc_start_main in "/lib/x86_64-linux-gnu/libc.so.6"
#11  ? at ??:?
Floating point exception (core dumped)


Divide diameter by 10


::

    $ constant/phaseProperties

        particles
        {
            diameterModel constant;
            constantCoeffs
            {
                d               27.5e-6;
            }

            alphaMax 0.62;
        }


Still the same - revert. Possible problem with the setFields - wrong dimensions


Try blockMesh
~~~~~~~~~~~~~

::

	$ rm -rf processor* 0.* [1-9]* postProcessing *.foam

::

	$ cp 0/alpha.air.org 0/alpha.air
	$ cp 0/alpha.particles.org 0/alpha.particles


::

	$ setFields

::

	$ twoPhaseEulerFoam > log &

::

	$ decomposePar

::

	$ mpirun -np 8 twoPhaseEulerFoam -parallel > log &


::

	$ kill -9 PID

::

	$ pkill twoPhaseEulerFo


Run for only 0.5 seconds

::

	$ nano system/controlDict

		endTime         0.5;


::

	$ reconstructPar


Modify blockMesh
~~~~~~~~~~~~~~~~


--> FOAM FATAL ERROR: 
face 0 in patch 1 does not have neighbour cell face: 4(11 15 14 10)
























Tutorial on Non-Newtonain fluid
===============================

fluentMeshToFoam needs a system directory and a controlDict:



::

	$ two
	$ mkdir system
	$ nano system/controlDict

		application     nonNewtonianIcoFoam;

		startFrom       startTime;

		startTime       0;

		stopAt          endTime;

		endTime         5;

		deltaT          0.0025;

		writeControl    runTime;

		writeInterval   0.05;

		purgeWrite      0;

		writeFormat     ascii;

		writePrecision  6;

		writeCompression off;

		timeFormat      general;

		timePrecision   6;

		runTimeModifiable true;

Writes the output to constant/polyMesh

::

	$ fluentMeshToFoam fluent.msh

Create a folder for boundary conditions:

::

	$ mkdir 0
	$ nano 0/p

		FoamFile
		{
			version     2.0;
			format      ascii;
			class       volScalarField;
			object      p;
		}

		dimensions      [0 2 -2 0 0 0 0];

		internalField   uniform 0;

		boundaryField
		{
			INLET
			{
				type            zeroGradient;
			}

			OUTLET
			{
				type            fixedValue;
				value           uniform 0;
			}


			FIXED_WALLS
			{
				type            zeroGradient;
			}
		

			FRONT_AND_BACK
			{
				type            empty;
			}

		}

	$ nano 0/U

		FoamFile
		{
			version     2.0;
			format      ascii;
			class       volVectorField;
			object      U;
		}

		dimensions      [0 1 -1 0 0 0 0];

		internalField   uniform (0 0 0);

		boundaryField
		{
			INLET
			{
				type            fixedValue;
				value           uniform (1 0 0);
			}

			OUTLET
			{
				type            zeroGradient;
			}


			FIXED_WALLS
			{
				type            fixedValue;
				value           uniform (0 0 0);
			}

		
			FRONT_AND_BACK
			{
				type            empty;
			}

		}

Copy BirdCarreau into directory:

::

	$ cd $WM_PROJECT_DIR
	$ cp -r --parents src/transportModels/incompressible/viscosityModels/BirdCarreau/ $WM_PROJECT_USER_DIR/
	$ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
	$ mv BirdCarreau Casson
	$ cd Casson
	$ mv BirdCarreau.C Casson.C
	$ mv BirdCarreau.H Casson.H
	$ sed -i s/BirdCarreau/Casson/g Casson.C
	$ sed –i s/BirdCarreau/Casson/g Casson.H			(may not work, so find and replace in file) 



Create make files:

::

	$ mkdir Make
	$ nano Make/files
		
		Casson.C

		LIB = $(FOAM_USER_LIBBIN)/libCasson

	$ nano Make/options

		EXE_INC = \
			-I$(LIB_SRC)/transportModels/incompressible/lnInclude \
			-I$(LIB_SRC)/finiteVolume/lnInclude
		LIB_LIBS = \
			-lfiniteVolume


Add formula for Casson plot:

::

	$ nano Casson.C

		// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

		Foam::tmp<Foam::volScalarField>
		Foam::viscosityModels::Casson::calcNu() const
		{
        // ADDED:
		return max
			(
				nuMin_, min( nuMax_,pow(pow(tau0_/max(strainRate(),dimensionedScalar("VSMALL", dimless/dimTime,VSMALL)),0.5)+pow(m_,0.5),scalar(2.0)))
			);

		}

Add new coefficients:

::

	// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

	Foam::viscosityModels::Casson::Casson
	(
		const word& name,
		const dictionary& viscosityProperties,
		const volVectorField& U,
		const surfaceScalarField& phi
	)
	:
		viscosityModel(name, viscosityProperties, U, phi),
		CassonCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
        m_(CassonCoeffs_.lookup("m")),            // ADDED
        tau0_(CassonCoeffs_.lookup("tau0")),      // ADDED
        nuMin_(CassonCoeffs_.lookup("nuMin")),    // ADDED
        nuMax_(CassonCoeffs_.lookup("nuMax")),    // ADDED
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

	bool Foam::viscosityModels::Casson::read
	(
		const dictionary& viscosityProperties
	)
	{
		viscosityModel::read(viscosityProperties);

		CassonCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

		CassonCoeffs_.lookup("nu0") >> nu0_;             // ADDED
		CassonCoeffs_.lookup("nuInf") >> nuInf_;         // ADDED
		CassonCoeffs_.lookup("k") >> k_;                 // ADDED
		CassonCoeffs_.lookup("n") >> n_;                 // ADDED


Add coefficients in header file:

::

	$ nano Casson.H
	
        dimensionedScalar m_;
        dimensionedScalar tau0_;
        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;

Compile in the main folder:

::

	$ wmake libso

Test the solver
---------------

In constant/transportProperties

::

	transportModel  Casson;

	CassonCoeffs
	{
		m             m [ 0 2 -1 0 0 0 0 ] 0.00414;
		tau0          tau0 [0 2 -2 0 0 0 0] 0.0038;
		nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
		nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
	}

In system/controlDict

::

	libs
	(
		"libCasson.so"
	);

::

	$ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
	$ foamCalc mag U

Sample the velocity profile:

::

	$ nano ./system/sampleDict

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
			object      sampleDict;
		}
		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

		interpolationScheme cellPoint;

		setFormat       raw;

		sets
		(
			data
			{
				type    uniform;
				axis    y;
				start   (10 0 0);
				end     (10 1 0);
				nPoints 10;
			}
		);

		fields          (magU);

		// ************************************************************************* //

Sample the data:

	$ sample

Run the non-Newtonian case

Remove the data:

::

	$ rm -rf 0.* [1-9]* postProcessing *.foam

In constant/transportProperties

::

	transportModel  Newtonian;

	nu              nu [ 0 2 -1 0 0 0 0 ] 1;

In constant/controlDict remove the library:

::

	libs
	(
		    "libCasson.so"
	);


::

	$ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
	$ foamCalc mag U

Sample the data:

	$ sample

Plot in Python:

::


	import numpy as np
	import matplotlib.pyplot as plt
	import math as mt

	data = np.genfromtxt('postProcessing/sets/5/data_magU.xy', delimiter=' ', skip_header=0, names=['x','magU'])
	data2 = np.genfromtxt('../non_newtonian_002_plates_test_casson/postProcessing/sets/5/data_magU.xy', delimiter=' ', skip_header=0, names=['x2','magU2'])
	x = data['x']
	magU = data['magU']

	x2 = data2['x2']
	magU2 = data2['magU2']

	ax1 = plt.gca()

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	line3=ax1.plot(data['x'], data['magU'], color='r', label='Newtonian')
	line4=ax1.plot(data2['x2'], data2['magU2'], color='g', label='Non-Newtonian')

	lines = line3+line4
	labels = [l.get_label() for l in lines]
	legend= ax1.legend(lines, labels, loc=0, fontsize='medium')

	ax1.set_xlabel(r'Distance, x \textit{(m)}')
	ax1.set_ylabel(r'Velocity Magnitude \textit{(m/s)}')
	ax1.set_xlim([0,1])
	ax1.set_ylim([0,2])
	plt.savefig('compare_velocity_profiles')
	plt.show()


Increase the mesh density near the wall - Newtonian
---------------------------------------------------

::

    $ cp -rf non_newtonian_003_plates_test_newtonian non_newtonian_004_plates_test_newtonian_dense
    $ cd non_newtonian_004_plates_test_newtonian_dense

::

	$ rm -rf 0.* [1-9]* postProcessing *.foam *.png

Increase mesh density

::

    $ nano constant/polyMesh/blockMeshDict
    
        blocks
        (
            hex (0 1 2 3 4 5 6 7) porosity (120 20 1)
            simpleGrading
            (
                1       // x-direction expansion
                (
                        (0.2 0.3 4)
                        (0.6 0.4 1)
                        (0.2 0.3 0.25)
                )
                1       // z-direction expansion
            )
        );

Re-run the nonNewtonian case:

::

	$ nonNewtonianIcoFoam
	
Calculate the velocity magnitude:

::
	
	$ foamCalc mag U



Sample the velocity profile:

::

	$ nano ./system/sampleDict

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
			object      sampleDict;
		}
		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

		interpolationScheme cellPoint;

		setFormat       raw;

		sets
		(
			data
			{
				type    uniform;
				axis    y;
				start   (10 0 0);
				end     (10 1 0);
				nPoints 1000;
			}
		);

		fields          (magU);

		// ************************************************************************* //

Sample the data:

	$ sample

Plot in Python:

::

    $ python plot.py
    
    
Increase the mesh density near the wall - Casson 
------------------------------------------------
    
::

    $ cp -rf non_newtonian_002_plates_test_casson non_newtonian_005_plates_test_casson_dense
    $ cp -rf non_newtonian_004_plates_test_newtonian_dense/constant/polyMesh/blockMeshDict non_newtonian_005_plates_test_casson_dense/constant/polyMesh/blockMeshDict
    $ cd non_newtonian_005_plates_test_casson_dense

::

	$ rm -rf 0.* [1-9]* postProcessing *.foam *.png
    
Re-run the nonNewtonian case:

::

	$ nonNewtonianIcoFoam
    
Sample the velocity profile:

::

	$ nano ./system/sampleDict

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
			object      sampleDict;
		}
		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

		interpolationScheme cellPoint;

		setFormat       raw;

		sets
		(
			data
			{
				type    uniform;
				axis    y;
				start   (10 0 0);
				end     (10 1 0);
				nPoints 1000;
			}
		);

		fields          (magU);

		// ************************************************************************* //
    
Sample the data:

::

    $ sample    
    
Change the python script:

::

    $ cd ../non_newtonian_004_plates_test_newtonian_dense    
    $ nano plot.py
        data2 = np.genfromtxt('../non_newtonian_005_plates_test_casson_dense/postProcessing/sets/5/data_magU.xy', delimiter=' ', skip_hea$
    $ python plot.py
    




Additional Modifications to Non-Newtonian fluid model
=====================================================




Add da Cruz non-Newtonian model
-------------------------------

* Add pressure as volScalarField (like temperature dependent viscosity)
* Set constant effective friction and nu = mu p / gamma
* Set linear effective friction

ConstantFriction

::

	$ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
	$ cp -rf Casson ConstantFriction
	$ cd ConstantFriction
	
Remove everything except .C and .H files:

::
	
	$ rm *.dep *.C.save
	$ mv Casson.C ConstantFriction.C
	$ mv Casson.H ConstantFriction.H
	$ sed -i s/Casson/ConstantFriction/g ConstantFriction.C
	$ sed -i s/Casson/ConstantFriction/g ConstantFriction.H			

Create make files:

::

	$ nano Make/files
		
		ConstantFriction.C

		LIB = $(FOAM_USER_LIBBIN)/libConstantFriction

	$ nano Make/options

		EXE_INC = \
			-I$(LIB_SRC)/transportModels/incompressible/lnInclude \
			-I$(LIB_SRC)/finiteVolume/lnInclude
		LIB_LIBS = \
			-lfiniteVolume

Add formula for ConstantFriction plot:

::

	$ nano ConstantFriction.C

        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        Foam::tmp<Foam::volScalarField>
        Foam::viscosityModels::ConstantFriction::calcNu() const
        {   
            const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");
            
            return max
            (
                nuMin_, 
                min
                (
                    nuMax_, 
                    (muStar_*p)/
                    max
                    (
                        strainRate(),
                        dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
                    )
                )
            );
        }


Add new coefficients:

::

	$ nano ConstantFriction.C
	
	    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        Foam::viscosityModels::ConstantFriction::ConstantFriction
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        )
        :
            viscosityModel(name, viscosityProperties, U, phi),
            ConstantFrictionCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
            muStar_(ConstantFrictionCoeffs_.lookup("muStar")),
            nuMin_(ConstantFrictionCoeffs_.lookup("nuMin")),
            nuMax_(ConstantFrictionCoeffs_.lookup("nuMax")),
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

::

	$ nano ConstantFriction.C
	
        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        bool Foam::viscosityModels::ConstantFriction::read
        (
            const dictionary& viscosityProperties
        )
        {
            viscosityModel::read(viscosityProperties);

            ConstantFrictionCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

            ConstantFrictionCoeffs_.lookup("muStar") >> muStar_;
            ConstantFrictionCoeffs_.lookup("nuMin") >> nuMin_;
            ConstantFrictionCoeffs_.lookup("nuMax") >> nuMax_;


            return true;
        }


        // ************************************************************************* //

Add coefficients in header file:

::

    $ nano ConstantFriction.H
	
        dimensionedScalar muStar_;
        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;


Compile in the main folder (creates .dep):

::

    $ wmake libso

If a mistake is made, clean using

::

    $ wclean libso

Test the solver
---------------

::

    $ cp -rf non_newtonian_005_plates_test_casson_dense non_newtonian_006_plates_test_constant_friction_dense
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

In constant/transportProperties

::

    $ nano constant/transportProperties
    
        transportModel  ConstantFriction;

        ConstantFrictionCoeffs
        {
                muStar        muStar [ 0 0 0 0 0 0 0 ] 0.11;
                nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
                nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
        }

In system/controlDict

::

    $ nano system/controlDict
    
        libs
        (
                "libConstantFriction.so"
        );
        
::

	$ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
	$ foamCalc mag U

Sample the data:

::

	$ sample    

Copy python

::

    $ cp ../non_newtonian_004_plates_test_newtonian_dense/plot.py .
	
Change the python script:

::

    $ cd ../non_newtonian_004_plates_test_newtonian_dense    
    $ nano plot.py
        data2 = np.genfromtxt('../non_newtonian_005_plates_test_casson_dense/postProcessing/sets/5/data_magU.xy', delimiter=' ', skip_hea$
    $ python plot.py  
        

Add linear daCruz model - LinearFriction
----------------------------------------

LinearFriction

::

	$ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
	$ cp -rf ConstantFriction LinearFriction
	$ cd LinearFriction
	
Remove everything except .C and .H files:

::
	
	$ wclean libso
	$ mv ConstantFriction.C LinearFriction.C
	$ mv ConstantFriction.H LinearFriction.H
	$ sed -i s/ConstantFriction/LinearFriction/g LinearFriction.C
	$ sed -i s/ConstantFriction/LinearFriction/g LinearFriction.H


Create make files:

::

	$ nano Make/files
		
		LinearFriction.C

		LIB = $(FOAM_USER_LIBBIN)/libLinearFriction

	$ nano Make/options

		EXE_INC = \
			-I$(LIB_SRC)/transportModels/incompressible/lnInclude \
			-I$(LIB_SRC)/finiteVolume/lnInclude
		LIB_LIBS = \
			-lfiniteVolume	
	

Add Linear Friction model:		

	
::

    $ nano LinearFriction.C	
	
	
        / * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        Foam::tmp<Foam::volScalarField>
        Foam::viscosityModels::LinearFriction::calcNu() const
        {
            const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");
            const volScalarField& muStar = muStar1_ + muStar2_ * (strainRate() / 
                                pow
                                (
                                    max
                                    (
                                        p,
                                        dimensionedScalar("VSMALL", p.dimensions(), VSMALL)
                                    ),
                                    0.5
                                )
                            );

        
            return max
            (
                nuMin_,
                min
                (
                    nuMax_,
                    (muStar*p)/
                    max
                    (
                        strainRate(),
                        dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
                    )
                )
            );

        }


Add new coefficients:

::

	$ nano LinearFriction.C
	
	    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        Foam::viscosityModels::LinearFriction::LinearFriction
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        )
        :
            viscosityModel(name, viscosityProperties, U, phi),
            LinearFrictionCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
            muStar1_(LinearFrictionCoeffs_.lookup("muStar1")),
            muStar2_(LinearFrictionCoeffs_.lookup("muStar2")),
            nuMin_(LinearFrictionCoeffs_.lookup("nuMin")),
            nuMax_(LinearFrictionCoeffs_.lookup("nuMax")),
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


::

	$ nano ConstantFriction.C
	
        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        bool Foam::viscosityModels::ConstantFriction::read
        (
            const dictionary& viscosityProperties
        )
        {
            viscosityModel::read(viscosityProperties);

            ConstantFrictionCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

            ConstantFrictionCoeffs_.lookup("muStar1") >> muStar1_;
            ConstantFrictionCoeffs_.lookup("muStar2") >> muStar2_; 
            ConstantFrictionCoeffs_.lookup("nuMin") >> nuMin_;
            ConstantFrictionCoeffs_.lookup("nuMax") >> nuMax_;


            return true;
        }


        // ************************************************************************* //

        
 Add coefficients in header file:

::

    $ nano LinearFriction.H
	
        dimensionedScalar muStar1_;
        dimensionedScalar muStar2_; 
        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;
       
Compile in the main folder (creates .dep):

::

    $ wmake libso

If a mistake is made, clean using

::

    $ wclean libso        
        
        
Test the solver
---------------

::

    $ cp -rf non_newtonian_006_plates_test_constant_friction_dense non_newtonian_007_plates_test_linear_friction_dense
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

In constant/transportProperties

::

    $ nano constant/transportProperties
    
        transportModel  LinearFriction;

        LinearFrictionCoeffs
        {
            muStar1       muStar1 [ 0 0 0 0 0 0 0 ] 0.11;
            muStar2       muStar2 [ 0 1 0 0 0 0 0 ] 1.62;
            nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
            nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
        }

In system/controlDict

::

    $ nano system/controlDict
    
        libs
        (
                "libLinearFriction.so"
        );
        
::

	$ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample    

Copy python

::

    $ cp ../non_newtonian_004_plates_test_newtonian_dense/plot.py .
	
Change the python script:

::

    $ cd ../non_newtonian_004_plates_test_newtonian_dense    
    $ nano plot.py
        data2 = np.genfromtxt('../non_newtonian_005_plates_test_casson_dense/postProcessing/sets/5/data_magU.xy', delimiter=' ', skip_hea$
    $ python plot.py          
        
        
::  

    Selecting incompressible transport model LinearFriction
    new cannot satisfy memory request.
    This does not necessarily mean you have run out of virtual memory.
    It could be due to a stack violation caused by e.g. bad use of pointers or an out of date shared library
    Aborted (core dumped)
        
        
::

    $ nano LinearFriction.C	
	
	
        / * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        Foam::tmp<Foam::volScalarField>
        Foam::viscosityModels::LinearFriction::calcNu() const
        {
            const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");
        
            return max
            (
                nuMin_,
                min
                (
                    nuMax_,
                    (( muStar1_ + muStar2_ * (strainRate() / pow (max(p,dimensionedScalar("VSMALL", p.dimensions(), VSMALL)),0.5)))*p) /
                    max
                    (
                        strainRate(),
                        dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
                    )
                )
            );

        }

        
Add non-linear Jop model - NonLinearFriction
--------------------------------------------        
        
NonLinearFriction

::

	$ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
	$ cp -rf LinearFriction NonLinearFriction
	$ cd NonLinearFriction
	
Remove everything except .C and .H files:

::

	$ two
	$ wclean libso
	$ mv LinearFriction.C NonLinearFriction.C
	$ mv LinearFriction.H NonLinearFriction.H
	$ sed -i s/LinearFriction/NonLinearFriction/g NonLinearFriction.C
	$ sed -i s/LinearFriction/NonLinearFriction/g NonLinearFriction.H
        
Create make files:

::

	$ nano Make/files
		
		NonLinearFriction.C

		LIB = $(FOAM_USER_LIBBIN)/libNonLinearFriction

	$ nano Make/options

		EXE_INC = \
			-I$(LIB_SRC)/transportModels/incompressible/lnInclude \
			-I$(LIB_SRC)/finiteVolume/lnInclude
		LIB_LIBS = \
			-lfiniteVolume	
	        
Add NonLinear Friction model:		

	
::

    $ nano NonLinearFriction.C	
	
	
        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

		Foam::tmp<Foam::volScalarField>
		Foam::viscosityModels::NonLinearFriction::calcNu() const
		{
			const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");

			return max
			(
				nuMin_,
				min
				(
				    nuMax_,
				    (( muStar1_ + ( muStar2_ - muStar1_ ) / ( iZero_ / (max(strainRate(),dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)) / pow (max(p,dimensionedScalar("VSMALL", p.dimensions(), VSMALL)),0.5)) + 1 ) )*p)/
				    max
				    (
				        strainRate(),
				        dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
				    )
				)
			);

		}
        
        
::

	$ nano LinearFriction.C
	
            // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

            Foam::viscosityModels::NonLinearFriction::NonLinearFriction
            (
                const word& name,
                const dictionary& viscosityProperties,
                const volVectorField& U,
                const surfaceScalarField& phi
            )
            :
                viscosityModel(name, viscosityProperties, U, phi),
                NonLinearFrictionCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
                iZero_(NonLinearFrictionCoeffs_.lookup("iZero")),
                muStar1_(NonLinearFrictionCoeffs_.lookup("muStar1")),
                muStar2_(NonLinearFrictionCoeffs_.lookup("muStar2")),
                nuMin_(NonLinearFrictionCoeffs_.lookup("nuMin")),
                nuMax_(NonLinearFrictionCoeffs_.lookup("nuMax")),
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

::            
            
    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    bool Foam::viscosityModels::NonLinearFriction::read
    (
        const dictionary& viscosityProperties
    )
    {
        viscosityModel::read(viscosityProperties);

        NonLinearFrictionCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");
        
        NonLinearFrictionCoeffs_.lookup("iZero") >> iZero_;
        NonLinearFrictionCoeffs_.lookup("muStar1") >> muStar1_;
        NonLinearFrictionCoeffs_.lookup("muStar2") >> muStar2_;
        NonLinearFrictionCoeffs_.lookup("nuMin") >> nuMin_;
        NonLinearFrictionCoeffs_.lookup("nuMax") >> nuMax_;


        return true;
    }            
            
Add coefficients in header file:

::

    $ nano NonLinearFriction.H
    
	dimensionedScalar iZero_;
        dimensionedScalar muStar1_;
        dimensionedScalar muStar2_; 
        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;            
            
 Compile in the main folder (creates .dep):

::

    $ wmake libso

If a mistake is made, clean using

::

    $ wclean libso             

    
Test the solver
---------------

::

    $ cp -rf non_newtonian_007_plates_test_linear_friction_dense non_newtonian_008_plates_test_non_linear_friction_dense
    $ cd non_newtonian_008_plates_test_non_linear_friction_dense
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

In constant/transportProperties

::

    $ nano constant/transportProperties
    
        transportModel  NonLinearFriction;

        NonLinearFrictionCoeffs
        {
            iZero         iZero [ 0 -1 0 0 0 0 0 ] 0.279;
            muStar1       muStar1 [ 0 0 0 0 0 0 0 ] 0.382;
            muStar2       muStar2 [ 0 0 0 0 0 0 0 ] 0.643;
            nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
            nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
        }

In system/controlDict

::

    $ nano system/controlDict
    
        libs
        (
                "libNonLinearFriction.so"
        );
        
::

	$ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample    
    
::

    Floating point exception (core dumped)
    
Made sure everything was positive

Now there is some distrubance - increase the length of the pipe to see if it is fully developed

::

    $ nano constant/polyMesh/blockMeshDict

        // The vertices can be defined using variables Lx, Ly, Lz. It saves time to modify the size of the domain

        Lx 20;
        Ly 1;
        Lz 0.1;

::

    $ nano system/sampleDict

        start   (20 0 0);
        end     (20 1 0);


        
::

    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

Create the new mesh:

::

    $ blockMesh    
    
    
::

    $ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample      

    
The inertial number is possibly too small. so that the regime being tested is not the dense regime.
So export the local strain rate to see what the inertial number is

::

    $ nano NonLinearFriction.C	
    
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
    ),
    strainRate_
    (
        IOobject
        (
            "strainRate",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        strainRate()
    )
    
::

    $ nano NonLinearFriction.H

        volScalarField strainRate_;
    
::

    $ cp -rf  non_newtonian_008_plates_test_non_linear_friction_dense non_newtonian_009_plates_test_non_linear_friction_dense_strain
    $ cd non_newtonian_009_plates_test_non_linear_friction_dense_strain
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png  
    
    
::

    $ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample  
    
    
Strain rate can't be zero, so output strain rate from a Functions

::

    $ nano NonLinearFriction.C

        Foam::tmp<Foam::volScalarField>
        Foam::viscosityModels::NonLinearFriction::calcStrainRate() const
        {

            return max
            (
                strainRate(),
                dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)

            );

        }
        
::

    $ nano NonLinearFriction.H
        
        //- Calculate and return the strain rate
        tmp<volScalarField> calcStrainRate() const;  
        
        
Strain rate is zero for most of the channel

So what is the friction?

::

    Foam::tmp<Foam::volScalarField>
    Foam::viscosityModels::NonLinearFriction::calcMu() const
    {

        return
        ( muStar1_ + ( muStar2_ - muStar1_ ) / ( iZero_ / (max(strainRate(),dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)) / pow (max(p,dimensionedScalar("VSMALL", p.dimensions(), VSMALL)),0.5)) + 1 ) );

    }

    mu_
    (
        IOobject
        (
            "mu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcMu()
    )
    
Viscosity is out of range because pressure is extremely high.

set pressure = 100 at inlet

::

    $ cp -rf  non_newtonian_009_plates_test_non_linear_friction_dense_strain non_newtonian_010_plates_non_linear_pressure
    $ cd non_newtonian_010_plates_non_linear_pressure
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png
    
    
    
    
::

	$ mkdir 0
	$ nano 0/p

            FoamFile
            {
                version     2.0;
                format      ascii;
                class       volScalarField;
                object      p;
            }
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

            dimensions      [0 2 -2 0 0 0 0];

            internalField   uniform 0;

            boundaryField
            {
                inlet
                {
                    type            fixedValue;
                    value           uniform 101.325;
                }

                outlet
                {
                    type            fixedValue;
                    value           uniform 0;
                }


                walls
                {
                    type            zeroGradient;
                }


                frontAndBack
                {
                    type            empty;
                }

            }


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
                    type            zeroGradient;
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
    
 ::

    $ nonNewtonianIcoFoam

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample  
       
Reduce pressure at inlet:       
              
::

    inlet
    {
        type            fixedValue;
        value           uniform 10;
    }    
    
Output in solver  not constructor for updated values
    
$ cp -r $FOAM_APP/solvers/incompressible/nonNewtonianIcoFoam nonNewtonianIcoFoamGranular

Once the laplacianFoam solver has been copied into the user directory, we rename the main file and edit the Make/files    
    
::    
    
    $ cd nonNewtonianIcoFoamGranular
    $ mv nonNewtonianIcoFoam.C nonNewtonianIcoFoamGranular.C
    $ wclean
    $ gedit Make/files

Change files    
    
::

    nonNewtonianIcoFoamGranular.C

    EXE = $(FOAM_USER_APPBIN)/nonNewtonianIcoFoamGranular

We can now clean the previous compilation with wclean and compile this new program with wmake.

::

    $ wclean
    $ wmake
    
Add this to NonLinearFriction.H

::

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
            strainRate_ = calcStrainRate();
            mu_ = calcMu();
            Info << "max(nu): " << max(calcNu()).value() << " min(nu): " << min(calcNu()).value() << endl; 
            Info << "max(sr): " << max(calcStrainRate()).value() << " min(sr): " << min(calcStrainRate()).value() << endl;
            Info << "max(mu): " << max(calcMu()).value() << " min(mu): " << min(calcMu()).value() << endl;
        }

        
        
        
It's unclear why the maximum value of calcstrainRate() as output from correct() is different from the value in paraview

before PISO loop:

max(calcStrainRate()).value() = 0.103329

after PISO loop:

max(calcStrainRate()).value()  = 0.0752414

Paraview: max strainrate (volumes) = 0.0629
                         (points) = 0.0833

May depend on the way the maximum is computed, using points or volumes     

Choose to define inertia number based on global values
    
Repeat with Newtonian ,  linear and Casson (Poiseulle)




non_newtonian_010_plates_non_linear_pressure        
--------------------------------------------

In system/controlDict

::

    $ nano system/controlDict
    
        application nonNewtonianIcoFoamGranular
        
::

    $ nonNewtonianIcoFoamGranular

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample
    
    
non_newtonian_011_plates_linear_pressure     
----------------------------------------

::

    $ cp -rf non_newtonian_010_plates_non_linear_pressure non_newtonian_011_plates_linear_pressure
    $ cd non_newtonian_011_plates_linear_pressure
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

In system/controlDict

::

    $ nano system/controlDict

    
        application nonNewtonianIcoFoamGranular
        
        libs
        (
                "libLinearFriction.so"
        );    
    
::

    $ nano constant/transportProperties
    
    transportModel  LinearFriction;
    
::

    $ nonNewtonianIcoFoamGranular

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample    

non_newtonian_012_plates_casson_pressure      
----------------------------------------

::

    $ cp -rf non_newtonian_011_plates_linear_pressure non_newtonian_012_plates_casson_pressure
    $ cd non_newtonian_012_plates_casson_pressure
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

In system/controlDict

::

    $ nano system/controlDict
    
        libs
        (
                "libCasson.so"
        );    
    
::

    $ nano constant/transportProperties
    
    transportModel  Casson;    

::

    $ nonNewtonianIcoFoam    
    

non_newtonian_013_plates_newtonian_pressure      
-------------------------------------------

::

    $ cp -rf non_newtonian_012_plates_casson_pressure non_newtonian_013_plates_newtonian_pressure
    $ cd non_newtonian_013_plates_newtonian_pressure
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

::

    $ nano system/controlDict
    
      //  libs
      //  (
      //         "libCasson.so"
      //  );    
    
::

    $ nano constant/transportProperties
    
    transportModel  Newtonian;    
    
::

    $ nonNewtonianIcoFoam 


non_newtonian_014_plates_linear_new_coefficients     
------------------------------------------------

Use coefficients from daCruz 2008

::

    $ cp -rf non_newtonian_011_plates_linear_pressure non_newtonian_014_plates_linear_new_coefficients
    $ cd non_newtonian_014_plates_linear_new_coefficients
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

In system/controlDict

::

    $ nano system/controlDict

    
        LinearFrictionCoeffs
        {
            muStar1       muStar1 [ 0 0 0 0 0 0 0 ] 0.22;
            muStar2       muStar2 [ 0 1 0 0 0 0 0 ] 1.00;
            nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
            nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
        }   
    
::

    $ nonNewtonianIcoFoamGranular

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample    
    
    
mu is becoming very large = 7.81189e+148, this is no -physical

Introduce satruation for mu

New Viscosity Model - LinearFrictionSaturated
---------------------------------------------

LinearFrictionSaturated

::
        
        $ two
        $ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
        $ cp -rf LinearFriction LinearFrictionSaturated
        $ cd LinearFrictionSaturated
	
Remove everything except .C and .H files:

::
	
	$ wclean libso
	$ mv LinearFriction.C LinearFrictionSaturated.C
	$ mv LinearFriction.H LinearFrictionSaturated.H
	$ sed -i s/LinearFriction/LinearFrictionSaturated/g LinearFrictionSaturated.C
	$ sed -i s/LinearFriction/LinearFrictionSaturated/g LinearFrictionSaturated.H
	$ sed -i s/LinearFriction/LinearFrictionSaturated/g Make/files

Add Linear Friction model:		

::

        $ nano LinearFrictionSaturated.C
        
            Foam::tmp<Foam::volScalarField>
            Foam::viscosityModels::LinearFrictionSaturated::calcNu() const
            {
                const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");

                return max
                (
                    nuMin_,
                    min
                    (
                        nuMax_,
                        ((min
                        (
                            muMax_,
                            muStar1_ + muStar2_ * ((strainRate() * diameter_) / pow (max(p,dimensionedScalar("VSMALL", p.dimensions(), VSMALL)),0.5))
                        ))*p)/
                        max
                        (
                            strainRate(),
                            dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
                        )
                    )
                );

            }
        
        
            Foam::tmp<Foam::volScalarField>
            Foam::viscosityModels::LinearFrictionSaturated::calcMu() const
            {
                const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");

                return min 
                (
                    muMax_,
                    muStar1_ + muStar2_ * ((strainRate() * diameter_) / pow (max(p,dimensionedScalar("VSMALL", p.dimensions(), VSMALL)),0.5))
                );

            }

            
            muMax_(LinearFrictionSaturatedCoeffs_.lookup("muMax")),
            diameter_(LinearFrictionSaturatedCoeffs_.lookup("diameter")),
            
            LinearFrictionSaturatedCoeffs_.lookup("muMax") >> muMax_;
            LinearFrictionSaturatedCoeffs_.lookup("diameter") >> diameter_;
            
            
::

    $ nano LinearFrictionSaturated.H

        dimensionedScalar muMax_;
        dimensionedScalar diameter_;

Compile in the main folder (creates .dep):

::

    $ wmake libso

If a mistake is made, clean using

::

    $ wclean libso 


non_newtonian_015_plates_linear_saturated     
-----------------------------------------

::

    $ cd $FOAM_RUN
    $ cp -rf non_newtonian_014_plates_linear_new_coefficients non_newtonian_015_plates_linear_saturated
    $ cd non_newtonian_015_plates_linear_saturated
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

::

    $ nano constant/transportProperties

    
        LinearFrictionSaturatedCoeffs
        {
            muStar1       muStar1 [ 0 0 0 0 0 0 0 ] 0.22;
            muStar2       muStar2 [ 0 1 0 0 0 0 0 ] 1.00;
            muMax         muMax [ 0 0 0 0 0 0 0 ] 0.44;
            nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
            nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
            diameter      diameter [ 0 1 0 0 0 0 0 ] 1;
        }   
   
::

    $ nano system/controlDict

        libs
        (
                "libLinearFrictionSaturated.so"
        );   
   
::

    $ nonNewtonianIcoFoamGranular

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample 


Introduce diameter into Inertial Number in Non Linear Friction Model
--------------------------------------------------------------------

NonLinearFrictionDiameter

::
        
        $ two
        $ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
        $ cp -rf NonLinearFriction NonLinearFrictionDiameter
        $ cd NonLinearFrictionDiameter
	
Remove everything except .C and .H files:

::
	
	$ wclean libso
	$ mv NonLinearFriction.C NonLinearFrictionDiameter.C
	$ mv NonLinearFriction.H NonLinearFrictionDiameter.H
	$ sed -i s/NonLinearFriction/NonLinearFrictionDiameter/g NonLinearFrictionDiameter.C
	$ sed -i s/NonLinearFriction/NonLinearFrictionDiameter/g NonLinearFrictionDiameter.H
	$ sed -i s/NonLinearFriction/NonLinearFrictionDiameter/g Make/files

Add Linear Friction model:		

::

    Foam::tmp<Foam::volScalarField>
    Foam::viscosityModels::NonLinearFrictionDiameter::calcNu() const
    {
        const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");
        //volScalarField sR("StrainRate", strainRate());
        //Info << "Max volSR = " << max(sR) << "; min volSR = " << min(sR) << endl;
        //Info << "max(strainRate): " << max(strainRate()).value() << " min(strainRate): " << min(strainRate()).value() << endl; 
        
        return max
        (
            nuMin_,
            min
            (
                nuMax_,
                (( muStar1_ + ( muStar2_ - muStar1_ ) / ( iZero_ / ( ((max(strainRate(),dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)) * diameter_) / pow (max(p,dimensionedScalar("VSMALL", p.dimensions(), VSMALL)),0.5)) + 1 ) ) )*p)/
                max
                (
                    strainRate(),
                    dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
                )
            )
        );

    }
    
    
    Foam::tmp<Foam::volScalarField>
    Foam::viscosityModels::NonLinearFrictionDiameter::calcMu() const
    {
        const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");

        return
        ( muStar1_ + ( muStar2_ - muStar1_ ) / ( iZero_ / ( ((max(strainRate(),dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)) * diameter_) / pow (max(p,dimensionedScalar("VSMALL", p.dimensions(), VSMALL)),0.5)) + 1 ) ) );

    }
    
    diameter_(NonLinearFrictionDiameterCoeffs_.lookup("diameter")),
    NonLinearFrictionDiameterCoeffs_.lookup("diameter") >> diameter_;
    
    
::

    $ nano LinearFrictionSaturated.H
    
        dimensionedScalar diameter_;
        
Compile in the main folder (creates .dep):

::

    $ wmake libso

If a mistake is made, clean using

::

    $ wclean libso 
    
    
non_newtonian_016_plates_non_linear_diameter     
--------------------------------------------

::

    $ cd $FOAM_RUN
    $ cp -rf non_newtonian_015_plates_linear_saturated non_newtonian_016_plates_linear_diameter
    $ cd non_newtonian_016_plates_linear_diameter
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

::

    $ nano constant/transportProperties

    
        NonLinearFrictionDiameterCoeffs
        {
                iZero         iZero [ 0 0 0 0 0 0 0 ] 0.279;
                muStar1       muStar1 [ 0 0 0 0 0 0 0 ] 0.382;
                muStar2       muStar2 [ 0 0 0 0 0 0 0 ] 0.643;
                nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
                nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
                diameter      diameter [ 0 1 0 0 0 0 0 ] 1;
        }   
   
::

    $ nano system/controlDict

        libs
        (
                "libNonLinearFrictionDiameter.so"
        );   
   
::

    $ nonNewtonianIcoFoamGranular

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample     
    
Declare variables separately for easy readability
-------------------------------------------------


NonLinearFrictionDiameterDeclare

::
        
        $ two
        $ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
        $ cp -rf NonLinearFrictionDiameter NonLinearFrictionDiameterDeclare
        $ cd NonLinearFrictionDiameterDeclare

Remove everything except .C and .H files:

::
	
	$ wclean libso
	$ mv NonLinearFrictionDiameter.C NonLinearFrictionDiameterDeclare.C
	$ mv NonLinearFrictionDiameter.H NonLinearFrictionDiameterDeclare.H
	$ sed -i s/NonLinearFrictionDiameter/NonLinearFrictionDiameterDeclare/g NonLinearFrictionDiameterDeclare.C
	$ sed -i s/NonLinearFrictionDiameter/NonLinearFrictionDiameterDeclare/g NonLinearFrictionDiameterDeclare.H
	$ sed -i s/NonLinearFrictionDiameter/NonLinearFrictionDiameterDeclare/g Make/files

::

	// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

	Foam::tmp<Foam::volScalarField>
	Foam::viscosityModels::NonLinearFrictionDiameterDeclare::calcNu() const
	{
		const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");
		
		const volScalarField strainRate_
		( 
		    max(strainRate(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL))     // Zero strainRate -> Divide by zero in return statement
		);
		
		return max
		(
		    nuMin_, min ( nuMax_, (calcMu()*p)/strainRate_ )
		);

	}


	Foam::tmp<Foam::volScalarField>
	Foam::viscosityModels::NonLinearFrictionDiameterDeclare::calcMu() const
	{
		const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");
		
		// VSMALL = 1.0e-37
		
		const volScalarField strainRate_
		( 
		    max(strainRate(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL))     // Zero strainRate -> Zero iNumber_ -> Divide by zero in return statement             
		);
		
		const volScalarField pressure_
		( 
		    max(p, dimensionedScalar("VSMALL", p.dimensions(), VSMALL))                 // Zero pressure -> Sqrt of negative number in iNumber_
		);
		
		const volScalarField iNumber_
		( 
		    (strainRate_ * diameter_) / pow (pressure_, 0.5) 
		);
		
		return ( muStar1_ + ( muStar2_ - muStar1_ ) / ( (iZero_ / iNumber_) + 1 ) );

	}





Compile in the main folder (creates .dep):

::

    $ wmake libso

If a mistake is made, clean using

::

    $ wclean libso     
    
non_newtonian_017_plates_nonlinear_diameter_declare     
---------------------------------------------------

::

    $ cd $FOAM_RUN
    $ cp -rf non_newtonian_016_plates_nonlinear_diameter non_newtonian_017_plates_nonlinear_diameter_declare
    $ cd non_newtonian_017_plates_nonlinear_diameter_declare
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

::

    $ nano constant/transportProperties

    
        NonLinearFrictionDiameterDeclareCoeffs
        {
                iZero         iZero [ 0 0 0 0 0 0 0 ] 0.279;
                muStar1       muStar1 [ 0 0 0 0 0 0 0 ] 0.382;
                muStar2       muStar2 [ 0 0 0 0 0 0 0 ] 0.643;
                nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
                nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
                diameter      diameter [ 0 1 0 0 0 0 0 ] 1;
        }   
   
::

    $ nano system/controlDict

        libs
        (
                "libNonLinearFrictionDiameterDeclare.so"
        );   
   
::

    $ nonNewtonianIcoFoamGranular

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample 

Declare variables separately for easy readability
-------------------------------------------------


LinearFrictionDiameterDeclare

::
        
    $ two
    $ cd $WM_PROJECT_USER_DIR/src/transportModels/incompressible/viscosityModels
    $ cp -rf LinearFrictionSaturated LinearFrictionDiameterDeclare
    $ cd LinearFrictionDiameterDeclare

Remove everything except .C and .H files:

::
	
	$ mv LinearFrictionSaturated.C LinearFrictionDiameterDeclare.C
	$ mv LinearFrictionSaturated.H LinearFrictionDiameterDeclare.H
	$ sed -i s/LinearFrictionSaturated/LinearFrictionDiameterDeclare/g LinearFrictionDiameterDeclare.C
	$ sed -i s/LinearFrictionSaturated/LinearFrictionDiameterDeclare/g LinearFrictionDiameterDeclare.H
	$ sed -i s/LinearFrictionSaturated/LinearFrictionDiameterDeclare/g Make/files

::

	Foam::tmp<Foam::volScalarField>
	Foam::viscosityModels::LinearFrictionDiameterDeclare::calcNu() const
	{
		const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");

		const volScalarField strainRate_
		( 
		    max(strainRate(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL))     // Zero strainRate -> Divide by zero in return statement
		);
		
		return max
		(
		    nuMin_, min ( nuMax_, (calcMu()*p)/strainRate_ )
		);

	}


	Foam::tmp<Foam::volScalarField>
	Foam::viscosityModels::LinearFrictionDiameterDeclare::calcMu() const
	{
		const volScalarField& p = U_.mesh().lookupObject<volScalarField>("p");

		const volScalarField pressure_
		( 
		    max(p, dimensionedScalar("VSMALL", p.dimensions(), VSMALL))                 // Zero pressure -> Sqrt of negative number in iNumber_
		);
		
		const volScalarField iNumber_
		( 
		    (strainRate() * diameter_) / pow (pressure_, 0.5) 
		);
		
		
		return min 
		(
		    muMax_, ( muStar1_ + muStar2_ * iNumber_ )
		);

	}


non_newtonian_018_plates_linear_diameter_declare     
------------------------------------------------

::

    $ cd $FOAM_RUN
    $ cp -rf non_newtonian_017_plates_nonlinear_diameter non_newtonian_018_plates_linear_diameter_declare
    $ cd non_newtonian_018_plates_linear_diameter_declare
    $ rm -rf 0.* [1-9]* postProcessing *.foam *.png

::

    $ nano constant/transportProperties

    
		LinearFrictionDiameterDeclareCoeffs
		{
			muStar1       muStar1 [ 0 0 0 0 0 0 0 ] 0.22;
			muStar2       muStar2 [ 0 0 0 0 0 0 0 ] 1.00;
			muMax         muMax [ 0 0 0 0 0 0 0 ] 0.44;
			nuMin         nuMin [ 0 2 -1 0 0 0 0 ] 0.0001;
			nuMax         nuMax [ 0 2 -1 0 0 0 0 ] 100;
			diameter      diameter [ 0 1 0 0 0 0 0 ] 1.0;
		}
 
   
::

    $ nano system/controlDict

        libs
        (
                "libLinearFrictionDiameterDeclare.so"
        );   
   
::

    $ nonNewtonianIcoFoamGranular

Calculate the velocity magnitude:

::
	
    $ foamCalc mag U

Sample the data:

::

    $ sample 


twoPhaseEulerFoam Modifications
===============================

::

	$ cp -rf /home/apr207/OpenFOAM/OpenFOAM-2.4.0/tutorials/multiphase/twoPhaseEulerFoam/laminar .

Laminar air flow, turbulent particle flow

Altered geometry so that the width is only 8cm (from 15cm)

Also added alpha.air to initial conditions so that setFields will work correctly

::

	$ cd ~/OpenFOAM/apr207-2.4.0/run/eulerian_006_fluidisedBed_geometry/

::

	$ ./script_run

::

	$ gnuplot -persist plot_residuals_live

control + c will quit this window

Wall clock time = 315 seconds (5 min)


Copy twoPhaseEulerFoam to check if 2.2.2 is installed correctly
---------------------------------------------------------------

Copy the code from source

::

	cd $WM_PROJECT_USER_DIR/applications/solvers
	cp -rf $FOAM_APP/solvers/multiphase/twoPhaseEulerFoam/ .
	mv twoPhaseEulerFoam twoPhaseEulerFoamNew
	cd twoPhaseEulerFoamNew
	mv twoPhaseEulerFoam.C twoPhaseEulerFoamNew.C
	sed -i s/twoPhaseEulerFoam/twoPhaseEulerFoamNew/g twoPhaseEulerFoamNew.C
	sed -i s/twoPhaseEulerFoam/twoPhaseEulerFoamNew/g Make/files
	wclean
	wmake

	In file included from twoPhaseEulerFoamNew.C:64:0:
	/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude/readTimeControls.H: In function ‘int main(int, char**)’:
	/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude/readTimeControls.H:38:8: warning: unused variable ‘maxDeltaT’ [-Wunused-variable]
	 scalar maxDeltaT =


Run the case:

::

	run
	cp -rf $FOAM_TUTORIALS/multiphase/twoPhaseEulerFoam/bed2/ .
	cd bed2
	sed -i s/twoPhaseEulerFoam/twoPhaseEulerFoamNew/g system/controlDict
	blockMesh


Corrected 2.1.0 case and solver to 2.2.2
----------------------------------------

Run the case:

::

	run
	cd twoPhaseEulerSource
	blockMesh
	cp 0/alpha1.org 0/alpha1
	setFields
	decomposePar
	mpirun -np 8 twoPhaseEulerSource -parallel > log &



Adding New Boundary Condition to twoPhaseEulerFoam
--------------------------------------------------

Create new solver
~~~~~~~~~~~~~~~~~

::

	cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase
	cp -rf twoPhaseEulerFoamNew twoPhaseEulerSrivastava
	cd twoPhaseEulerSrivastava
	mv twoPhaseEulerFoamNew.C twoPhaseEulerSrivastava.C
	sed -i s/twoPhaseEulerFoamNew/twoPhaseEulerSrivastava/g twoPhaseEulerSrivastava.C
	sed -i s/twoPhaseEulerFoamNew/twoPhaseEulerSrivastava/g Make/files
	wclean
	wmake


Copy JohnsonJacksonParticleSlip from 2.3.x to 2.2.2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenFOAM-2.3.x/applications/solvers/multiphase/twoPhaseEulerFoam/phaseCompressibleTurbulenceModels/kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleSlip/

into folder 

~/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerSrivastava/kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleSlip/

Added .C file to kineticTheoryModels/Make

::

    derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.C

    
Correct Make and options files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

twoPhaseEulerSrivastava    
""""""""""""""""""""""" 

nano Make/files (altered) - location and name of executable has changed

The first line means the .C file to compile is located one directory above the Make folder.

The second line means put the executable in the directory $FOAM_USER_APPBIN

::

    twoPhaseEulerSrivastava.C

    EXE = $(FOAM_USER_APPBIN)/twoPhaseEulerSrivastava

https://cfd.direct/openfoam/user-guide/compiling-applications/    
    
nano Make/options (altered):

* EXE_INC is the path to the libraries
* EXE_LIBS is the library names (these have been altered). The library path must also be given, the user location.
* The actual library files to be linked must be specified using the -l option and removing the lib prefix and .so extension from the library file name, e.g.  libnew.so is included with the flag -lnew

Why is -L specified only once? Possibly because -L adds an additional path to the default path ($FOAM_LIBBIN)

::

    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
        -IturbulenceModel \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastava/kineticTheoryModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastava/interfacialModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastava/phaseModel/lnInclude \
        -Iaveraging

    EXE_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lEulerianInterfacialModelsSrivastava \
        -lfiniteVolume \
        -lmeshTools \
        -lincompressibleTransportModels \
        -lphaseModelSrivastava \
        -lkineticTheoryModelSrivastava

The result:

::

    + wmake
    SOURCE=twoPhaseEulerSrivastava.C ;  
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude 
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/transportModels/incompressible/lnInclude 
    -IturbulenceModel 
    -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerSrivastava/kineticTheoryModels/lnInclude 
    -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerSrivastava/interfacialModels/lnInclude 
    -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerSrivastava/phaseModel/lnInclude 
    -Iaveraging 

    -L/home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib 
    -lEulerianInterfacialModelsSrivastava 
    -lfiniteVolume 
    -lmeshTools 
    -lincompressibleTransportModels 
    -lphaseModelSrivastava 
    -lkineticTheoryModelSrivastava 

Interfacial Model    
"""""""""""""""""             

The location of the nine .C files to compile is the same relative to the Make folder.

nano interfacialModels/Make/files (altered):

* The location and name of the EulerianInterfacialModels library has changed

::

    dragModels/dragModel/dragModel.C
    dragModels/dragModel/newDragModel.C
    dragModels/Ergun/Ergun.C
    dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C
    dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C
    dragModels/SchillerNaumann/SchillerNaumann.C
    dragModels/Gibilaro/Gibilaro.C
    dragModels/WenYu/WenYu.C
    dragModels/SyamlalOBrien/SyamlalOBrien.C

    LIB = $(FOAM_USER_LIBBIN)/libEulerianInterfacialModelsSrivastava

nano interfacialModels/Make/options (altered):

* The name and location of the phaseModel library has changed
    
::

    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastava/phaseModel/lnInclude

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lphaseModelSrivastava  
    
    
The result (just the bits that correspond with the above) - library `lphaseModelSrivastava` can't be seen - for some unknown reason. 

::

    + wmake libso interfacialModels
    SOURCE=dragModels/dragModel/dragModel.C ;  
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude 
    -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerSrivastava/phaseModel/lnInclude 
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude 
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   
    SOURCE=dragModels/dragModel/newDragModel.C ;  
    ...   
    SOURCE=dragModels/Ergun/Ergun.C ;  
    ...
    SOURCE=dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C ; 
    ...
    SOURCE=dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C ; 
    ...
    SOURCE=dragModels/SchillerNaumann/SchillerNaumann.C ;  
    ...
    SOURCE=dragModels/WenYu/WenYu.C ;  
    ...
    SOURCE=dragModels/Gibilaro/Gibilaro.C ; 
    ...
    SOURCE=dragModels/SyamlalOBrien/SyamlalOBrien.C ;  
    ...
    '/home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib/libEulerianInterfacialModelsSrivastava.so' is up to date.

    
Kinetic Theory Model    
""""""""""""""""""""      

nano kineticTheoryModel/Make/files (altered):

* The location and name of the kineticTheoryModel library has changed
    
::

    kineticTheoryModel/kineticTheoryModel.C

    viscosityModel/viscosityModel/viscosityModel.C
    viscosityModel/viscosityModel/newViscosityModel.C
    viscosityModel/Gidaspow/GidaspowViscosity.C
    viscosityModel/Syamlal/SyamlalViscosity.C
    viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C
    viscosityModel/none/noneViscosity.C

    conductivityModel/conductivityModel/conductivityModel.C
    conductivityModel/conductivityModel/newConductivityModel.C
    conductivityModel/Gidaspow/GidaspowConductivity.C
    conductivityModel/Syamlal/SyamlalConductivity.C
    conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C

    radialModel/radialModel/radialModel.C
    radialModel/radialModel/newRadialModel.C
    radialModel/CarnahanStarling/CarnahanStarlingRadial.C
    radialModel/LunSavage/LunSavageRadial.C
    radialModel/SinclairJackson/SinclairJacksonRadial.C

    granularPressureModel/granularPressureModel/granularPressureModel.C
    granularPressureModel/granularPressureModel/newGranularPressureModel.C
    granularPressureModel/Lun/LunPressure.C
    granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C

    frictionalStressModel/frictionalStressModel/frictionalStressModel.C
    frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C
    frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C
    frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C
    
    derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.C
    
    LIB = $(FOAM_USER_LIBBIN)/libkineticTheoryModelSrivastava

nano kineticTheoryModel/Make/options (altered):

::
    
    EXE_INC = \
        -I$(LIB_SRC)/foam/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastava/phaseModel/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastava/interfacialModels/lnInclude   


The result (I miss out the repeating block and any extra bits):
        
::

    + wmake libso kineticTheoryModels

    SOURCE=kineticTheoryModel/kineticTheoryModel.C

    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude 
    -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerSrivastava/phaseModel/lnInclude 
    -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerSrivastava/interfacialModels/lnInclude 
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude 
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude

    SOURCE=viscosityModel/viscosityModel/viscosityModel.C 
    SOURCE=viscosityModel/viscosityModel/newViscosityModel.C
    SOURCE=viscosityModel/Gidaspow/GidaspowViscosity.C ;  
    SOURCE=viscosityModel/Syamlal/SyamlalViscosity.C ;  
    SOURCE=viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C ;  
    SOURCE=conductivityModel/conductivityModel/conductivityModel.C ;  
    SOURCE=viscosityModel/none/noneViscosity.C ;  
    SOURCE=conductivityModel/conductivityModel/newConductivityModel.C ;  
    SOURCE=conductivityModel/Gidaspow/GidaspowConductivity.C ;  
    SOURCE=conductivityModel/Syamlal/SyamlalConductivity.C ;  
    SOURCE=conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C ; 
    SOURCE=radialModel/radialModel/radialModel.C ;  
    SOURCE=radialModel/radialModel/newRadialModel.C ;  
    SOURCE=radialModel/CarnahanStarling/CarnahanStarlingRadial.C ;  
    SOURCE=radialModel/LunSavage/LunSavageRadial.C ;  
    SOURCE=radialModel/SinclairJackson/SinclairJacksonRadial.C ;  
    SOURCE=granularPressureModel/granularPressureModel/granularPressureModel.C ;  
    SOURCE=granularPressureModel/granularPressureModel/newGranularPressureModel.C ; 
    SOURCE=granularPressureModel/Lun/LunPressure.C ;  
    SOURCE=granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C ; 
    SOURCE=frictionalStressModel/frictionalStressModel/frictionalStressModel.C ;  
    SOURCE=frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C ; 
    SOURCE=frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C ;  
    SOURCE=frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C ; 
    SOURCE=derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.C ;  
    '/home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib/libkineticTheoryModelSrivastava.so' is up to date.       
        
        
        
Phase Model    
"""""""""""

nano phaseModel/Make/files (altered)

* The location and name of the library has changed

::

    LIB = $(FOAM_USER_LIBBIN)/libphaseModelSrivastava

nano phaseModel/Make/options (unaltered)

::
    
    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude

    LIB_LIBS = \
        -lincompressibleTransportModels    

The result (just the bits that correspond with the above) - library `lincompressibleTransportModels` can't be seen:        
        
::

    + wmake libso phaseModel
    SOURCE=phaseModel/phaseModel.C
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude 
    -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/transportModels/incompressible/lnInclude 
    '/home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib/libphaseModelSrivastava.so' is up to date.
    
    
Test if case works
~~~~~~~~~~~~~~~~~~

Run the case:

::

    run
    cp -rf $FOAM_TUTORIALS/multiphase/twoPhaseEulerFoam/bed2/ .
    cd bed2
    sed -i s/twoPhaseEulerFoam/twoPhaseEulerSrivastava/g system/controlDict
    blockMesh
    twoPhaseEulerSrivastava


Is the partial slip condition being applied correctly?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes, it is being applied correctly, although it ignores friction


    
Consider comparing JohnsonJackson with JohnsonJacksonSchaeffer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Look to twophaseeulersedfoam for example from 2.1.0

         
OpenFOAM version 5.0
====================

OpenFOAM version 5.0 run script
-------------------------------

::

    # Run script for OpenFOAM 5.0

    . $WM_PROJECT_DIR/bin/tools/RunFunctions	                                    # loads RunFunctions

    rm -rf log.*

    export NP=6					                                    # number of processors

    # Mesh
    runApplication blockMesh

    # Check the mesh
    runApplication checkMesh

    # Decompose
    runApplication decomposePar

    # Run
    runParallel $(getApplication)

    # Reconstruct
    runApplication reconstructPar

    rm -r processor*

Rename surfaces    
---------------

::

    sed -i s/frontAndBackPlanes/frontAndBack/g alpha.particles
    sed -i s/frontAndBackPlanes/frontAndBack/g alpha.particles.orig
    sed -i s/frontAndBackPlanes/frontAndBack/g alphat.particles
    sed -i s/frontAndBackPlanes/frontAndBack/g epsilon.air
    sed -i s/frontAndBackPlanes/frontAndBack/g k.air
    sed -i s/frontAndBackPlanes/frontAndBack/g nut.air
    sed -i s/frontAndBackPlanes/frontAndBack/g nut.particles
    sed -i s/frontAndBackPlanes/frontAndBack/g p
    sed -i s/frontAndBackPlanes/frontAndBack/g p_rgh
    sed -i s/frontAndBackPlanes/frontAndBack/g T.air
    sed -i s/frontAndBackPlanes/frontAndBack/g T.particles
    sed -i s/frontAndBackPlanes/frontAndBack/g Theta.particles
    sed -i s/frontAndBackPlanes/frontAndBack/g U.air
    sed -i s/frontAndBackPlanes/frontAndBack/g U.particles

See if case run    
    
Set fields dictionary
---------------------

::

    defaultFieldValues
    (
        volScalarFieldValue alpha.air 1
        volScalarFieldValue alpha.particles 0
    );

    regions
    (
        boxToCell
        {
            box ( -0.04 0.0 0.0 ) ( 0.04 0.5 0.01 );
            fieldValues
            (
                volScalarFieldValue alpha.air 0.40
                volScalarFieldValue alpha.particles 0.60
            );
        }
    );

Turn off energy equation 
------------------------

Turn off e:

::

    "(h|e).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0;
        minIter         0;  // APR changed from 1 to 0
        maxIter         0;  // APR added
    }

Make air density constant
-------------------------
    
Change air properties:

::

    thermoType
    {
        type            heRhoThermo;
        mixture         pureMixture;
        transport       const;
        thermo          hConst;
        equationOfState rhoConst; // APR changed from perfectGas to rhoConst
        specie          specie;
        energy          sensibleInternalEnergy;
    }

    mixture
    {
        specie
        {
            molWeight   28.9;
        }
        equationOfState             // APR added equationOfState
        {
            rho         1.2;
        }
        thermodynamics
        {
            Cp          1007;
            Hf          0;
        }
        transport
        {
            mu          1.84e-05;
            Pr          0.7;
        }
    }
            
Execution time = 1631.89 seconds = 27 minutes


Phase properties
----------------

residualAlpha is the residual phase-fraction for a given phase. It is used to stabilize the phase momentum as the phase-fraction -> 0

blending method is for blended drag models "none" is a valid selection

sigma keyword is surface tension

The drag models offer swarm correction of the drag force, since it is observed that swarms of bubbles behave different from single bubbles. none indicates there is no correction (unity multiplies the drag coefficient). "none" is a valid selection.

No virtual mass

There has to be a heat transfer model - there is no option "none". "RanzMarshall" is a valid selection

::

    particles
    {
        diameterModel constant;
        constantCoeffs
        {
            d               100e-6;
        }

        alphaMax        0.65;
        residualAlpha   1e-6;
    }   
    
    
    drag
    (
        (particles in air)
        {
            type            WenYu;
            residualRe      1e-3;
            swarmCorrection
            {
                type        none;
            }
        }
    ); 
    
    virtualMass
    (
        (solids in gas)
        {
            type            none;       // APR changed from constantCoefficient
            //Cvm             0;        // APR changed from 0.5 to 0
        }
    );
    

thermophysicalProperties.air
----------------------------

The thermophysicalProperties dictionary is read by any solver that uses the thermophysical model library

heRhoThermo: for solvers that construct rhoThermo, such as twoPhaseEulerFoam

thermophysical models without reactions is pureMixture, which represents a mixture with fixed composition. When pureMixture is specified, the thermophysical models coefficients are specified within a sub-dictionary called mixture.

* transport model: const assumes a constant dynamic viscosity (mu) and Prandtl number (Pr)
* thermo model: assumes a constant specific heat capacity (Cp) and a heat of fusion (Hf)
* There is currently only one option for the specie model which specifies the composition of each constituent
* sensibleInternalEnergy is the form of the energy equation

::
    
    thermoType
    {
        type            heRhoThermo;
        mixture         pureMixture;
        transport       const;
        thermo          hConst;
        equationOfState rhoConst; // APR changed from perfectGas to rhoConst
        specie          specie;
        energy          sensibleInternalEnergy;
    }

    mixture
    {
        transport
        {
            mu          1.84e-05;
            Pr          0.7;
        }
        thermodynamics
        {
            Cp          1007;
            Hf          0;
        }
        equationOfState             // APR added equationOfState
        {
            rho         1.2;
        }
        specie
        {
            molWeight   28.9;
        }
    }   
    
Create location for version 5 solver
------------------------------------

::

    mkdir 


Install OpenFOAM 5.0 on rodrigo
-------------------------------

Download packages:

::

    download from http://dl.openfoam.org/source/5-0
    download from http://dl.openfoam.org/third-party/5-0
    scp -r OpenFOAM-5.x-version-5.0.tar.gz apr207@emps-rodrigo:/home/links/apr207/OpenFOAM
    scp -r ThirdParty-5.x-version-5.0.tar.gz apr207@emps-rodrigo:/home/links/apr207/OpenFOAM
    tar -zxvf OpenFOAM-5.x-version-5.0.tar.gz
    tar -zxvf ThirdParty-5.x-version-5.0.tar.gz
    mv OpenFOAM-5.x-version-5.0 OpenFOAM-5.0
    mv ThirdParty-5.x-version-5.0 ThirdParty-5.0
    
Don't need to install any software for compliation    
    
Set environment variables:    
    
::

    cd $HOME
    nano .bashrc
        alias five="module load mpi; source $HOME/OpenFOAM/OpenFOAM-5.0/etc/bashrc"

Install Third Party Scotch and PT Scotch
        
 ::
 
    five
    cd $WM_THIRD_PARTY_DIR
    ./Allwmake
    
Install OpenFOAM

::

    cd ../OpenFOAM-5.0/
    ./Allwmake
    
If you see this: 

::
        
    five    
        bash: mpicc: command not found...

Then:

::

    module load mpi
    
    
Copy across OpenFOAM 5.0 twoPhaseEulerFoam solver
-------------------------------------------------
 
 
 Copy the code from source

::

    five
    cd $WM_PROJECT_USER_DIR
    mkdir -p applications/solvers/multiphase
    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase
    cp -rf $FOAM_APP/solvers/multiphase/twoPhaseEulerFoam/ .
    mv twoPhaseEulerFoam twoPhaseEulerSrivastavaFoam
    cd twoPhaseEulerSrivastavaFoam
    mv twoPhaseEulerFoam.C twoPhaseEulerSrivastavaFoam.C
    sed -i s/twoPhaseEulerFoam/twoPhaseEulerSrivastavaFoam/g twoPhaseEulerSrivastavaFoam.C
    sed -i s/twoPhaseEulerFoam/twoPhaseEulerSrivastavaFoam/g Make/files
    ./Allwclean
    ./Allwmake

    
    
Correct Make and options files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

twoPhaseEulerSrivastava    
""""""""""""""""""""""" 

nano Make/files (altered) - location and name of executable has changed

The first line means the .C file to compile is located one directory above the Make folder.

The second line means put the executable in the directory $FOAM_USER_APPBIN

::

    twoPhaseEulerSrivastavaFoam.C

    EXE = $(FOAM_USER_APPBIN)/twoPhaseEulerSrivastavaFoam

https://cfd.direct/openfoam/user-guide/compiling-applications/    
    
nano Make/options (altered):

* EXE_INC is the path to the libraries
* EXE_LIBS is the library names (these have been altered). The library path must also be given, the user location.
* The actual library files to be linked must be specified using the -l option and removing the lib prefix and .so extension from the library file name, e.g.  libnew.so is included with the flag -lnew

Why is -L specified only once? Possibly because -L adds an additional path to the default path ($FOAM_LIBBIN)

::

    EXE_INC = \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/sampling/lnInclude

    EXE_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie \
        -lturbulenceModels \
        -lcompressibleTurbulenceModels \
        -lphaseCompressibleTurbulenceModelsSrivastava \
        -lincompressibleTransportModels \
        -lcompressibleTwoPhaseSystemSrivastava \
        -lcompressibleEulerianInterfacialModelsSrivastava \
        -lfiniteVolume \
        -lfvOptions \
        -lmeshTools \
        -lsampling

The result:



        
        
interfacialModels    
"""""""""""""""""             

The location of the nine .C files to compile is the same relative to the Make folder.

nano interfacialModels/Make/files (altered):

* The location and name of the EulerianInterfacialModels library has changed

::

    LIB = $(FOAM_USER_LIBBIN)/libcompressibleEulerianInterfacialModelsSrivastava

nano interfacialModels/Make/options (altered):

* The name and location of the phaseModel library has changed
    
::

    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/transportModel \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lcompressibleTwoPhaseSystemSrivastava \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie

        
phaseCompressibleTurbulenceModels
"""""""""""""""""""""""""""""""""

nano phaseCompressibleTurbulenceModels/Make/files (altered):

::

    LIB = $(FOAM_USER_LIBBIN)/libphaseCompressibleTurbulenceModelsSrivastava
       
        
nano phaseCompressibleTurbulenceModels/Make/options (altered):        

::

    EXE_INC = \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude\
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/transportModel \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie \
        -lturbulenceModels \
        -lcompressibleTurbulenceModels \
        -lincompressibleTransportModels \
        -lcompressibleTwoPhaseSystemSrivastava \
        -lcompressibleEulerianInterfacialModelsSrivastava \
        -lfiniteVolume \
        -lfvOptions \
        -lmeshTools
        
        
twoPhaseSystem
""""""""""""""

nano twoPhaseSystem/Make/files

::

     LIB = $(FOAM_USER_LIBBIN)/libcompressibleTwoPhaseSystemSrivastava   
  
  
nano twoPhaseSystem/Make/options  
  
::

    EXE_INC = \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/sampling/lnInclude

    LIB_LIBS = \
        -lincompressibleTransportModels \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie

Result:        
        
::

    wmake twoPhaseSystem
    wmakeLnInclude: linking include files to ./lnInclude

    -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem 
    -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude 
    -I/opt/openfoam5/src/transportModels/compressible/lnInclude 
    -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude 
    -I/opt/openfoam5/src/transportModels/incompressible/lnInclude 
    -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude 
    -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude 
    -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude 
    -I/opt/openfoam5/src/finiteVolume/lnInclude 
    -I/opt/openfoam5/src/meshTools/lnInclude 
    -I/opt/openfoam5/src/sampling/lnInclude 
    
        

Extended result:

::

    Allwmake /home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam
    wmakeLnInclude: linking include files to interfacialModels/lnInclude
    wmake twoPhaseSystem
    wmakeLnInclude: linking include files to ./lnInclude
    Making dependency list for source file twoPhaseSystem.C
    Making dependency list for source file orderedPhasePair.C
    Making dependency list for source file phasePair.C
    Making dependency list for source file phasePairKey.C
    Making dependency list for source file hyperbolic.C
    Making dependency list for source file linear.C
    Making dependency list for source file noBlending.C
    Making dependency list for source file newBlendingMethod.C
    Making dependency list for source file blendingMethod.C
    Making dependency list for source file randomCoalescence.C
    Making dependency list for source file turbulentBreakUp.C
    Making dependency list for source file wakeEntrainmentCoalescence.C
    Making dependency list for source file IATEsource.C
    Making dependency list for source file IATE.C
    Making dependency list for source file isothermalDiameter.C
    Making dependency list for source file constantDiameter.C
    Making dependency list for source file newDiameterModel.C
    Making dependency list for source file diameterModel.C
    Making dependency list for source file phaseModel.C
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c phaseModel/phaseModel.C -o Make/linux64GccDPInt32Opt/phaseModel/phaseModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/diameterModel/diameterModel.C -o Make/linux64GccDPInt32Opt/diameterModels/diameterModel/diameterModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/diameterModel/newDiameterModel.C -o Make/linux64GccDPInt32Opt/diameterModels/diameterModel/newDiameterModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/constantDiameter/constantDiameter.C -o Make/linux64GccDPInt32Opt/diameterModels/constantDiameter/constantDiameter.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/isothermalDiameter/isothermalDiameter.C -o Make/linux64GccDPInt32Opt/diameterModels/isothermalDiameter/isothermalDiameter.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/IATE/IATE.C -o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATE.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/IATE/IATEsources/IATEsource/IATEsource.C -o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATEsources/IATEsource/IATEsource.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/IATE/IATEsources/wakeEntrainmentCoalescence/wakeEntrainmentCoalescence.C -o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATEsources/wakeEntrainmentCoalescence/wakeEntrainmentCoalescence.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/IATE/IATEsources/turbulentBreakUp/turbulentBreakUp.C -o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATEsources/turbulentBreakUp/turbulentBreakUp.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c diameterModels/IATE/IATEsources/randomCoalescence/randomCoalescence.C -o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATEsources/randomCoalescence/randomCoalescence.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c BlendedInterfacialModel/blendingMethods/blendingMethod/blendingMethod.C -o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/blendingMethod/blendingMethod.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c BlendedInterfacialModel/blendingMethods/blendingMethod/newBlendingMethod.C -o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/blendingMethod/newBlendingMethod.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c BlendedInterfacialModel/blendingMethods/noBlending/noBlending.C -o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/noBlending/noBlending.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c BlendedInterfacialModel/blendingMethods/linear/linear.C -o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/linear/linear.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c BlendedInterfacialModel/blendingMethods/hyperbolic/hyperbolic.C -o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/hyperbolic/hyperbolic.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c phasePair/phasePairKey/phasePairKey.C -o Make/linux64GccDPInt32Opt/phasePair/phasePairKey/phasePairKey.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c phasePair/phasePair/phasePair.C -o Make/linux64GccDPInt32Opt/phasePair/phasePair/phasePair.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c phasePair/orderedPhasePair/orderedPhasePair.C -o Make/linux64GccDPInt32Opt/phasePair/orderedPhasePair/orderedPhasePair.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c twoPhaseSystem.C -o Make/linux64GccDPInt32Opt/twoPhaseSystem.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -shared -Xlinker --add-needed -Xlinker --no-as-needed Make/linux64GccDPInt32Opt/phaseModel/phaseModel.o Make/linux64GccDPInt32Opt/diameterModels/diameterModel/diameterModel.o Make/linux64GccDPInt32Opt/diameterModels/diameterModel/newDiameterModel.o Make/linux64GccDPInt32Opt/diameterModels/constantDiameter/constantDiameter.o Make/linux64GccDPInt32Opt/diameterModels/isothermalDiameter/isothermalDiameter.o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATE.o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATEsources/IATEsource/IATEsource.o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATEsources/wakeEntrainmentCoalescence/wakeEntrainmentCoalescence.o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATEsources/turbulentBreakUp/turbulentBreakUp.o Make/linux64GccDPInt32Opt/diameterModels/IATE/IATEsources/randomCoalescence/randomCoalescence.o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/blendingMethod/blendingMethod.o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/blendingMethod/newBlendingMethod.o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/noBlending/noBlending.o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/linear/linear.o Make/linux64GccDPInt32Opt/BlendedInterfacialModel/blendingMethods/hyperbolic/hyperbolic.o Make/linux64GccDPInt32Opt/phasePair/phasePairKey/phasePairKey.o Make/linux64GccDPInt32Opt/phasePair/phasePair/phasePair.o Make/linux64GccDPInt32Opt/phasePair/orderedPhasePair/orderedPhasePair.o Make/linux64GccDPInt32Opt/twoPhaseSystem.o -L/opt/openfoam5/platforms/linux64GccDPInt32Opt/lib \
        -lincompressibleTransportModels -lcompressibleTransportModels -lfluidThermophysicalModels -lspecie  -o /home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libcompressibleTwoPhaseSystemSrivastava.so
    wmake interfacialModels
    Making dependency list for source file wallDependentModel.C
    Making dependency list for source file Wellek.C
    Making dependency list for source file VakhrushevEfremov.C
    Making dependency list for source file TomiyamaAspectRatio.C
    Making dependency list for source file constantAspectRatio.C
    Making dependency list for source file newAspectRatioModel.C
    Making dependency list for source file aspectRatioModel.C
    Making dependency list for source file LopezDeBertodano.C
    Making dependency list for source file Gosman.C
    Making dependency list for source file Burns.C
    Making dependency list for source file constantTurbulentDispersionCoefficient.C
    Making dependency list for source file noTurbulentDispersion.C
    Making dependency list for source file newTurbulentDispersionModel.C
    Making dependency list for source file turbulentDispersionModel.C
    Making dependency list for source file TomiyamaWallLubrication.C
    Making dependency list for source file Frank.C
    Making dependency list for source file Antal.C
    Making dependency list for source file noWallLubrication.C
    Making dependency list for source file newWallLubricationModel.C
    Making dependency list for source file wallLubricationModel.C
    Making dependency list for source file Lamb.C
    Making dependency list for source file constantVirtualMassCoefficient.C
    Making dependency list for source file noVirtualMass.C
    Making dependency list for source file newVirtualMassModel.C
    Making dependency list for source file virtualMassModel.C
    Making dependency list for source file sphericalHeatTransfer.C
    Making dependency list for source file RanzMarshall.C
    Making dependency list for source file newHeatTransferModel.C
    Making dependency list for source file heatTransferModel.C
    Making dependency list for source file TomiyamaLift.C
    Making dependency list for source file LegendreMagnaudet.C
    Making dependency list for source file Moraga.C
    Making dependency list for source file constantLiftCoefficient.C
    Making dependency list for source file noLift.C
    Making dependency list for source file newLiftModel.C
    Making dependency list for source file liftModel.C
    Making dependency list for source file TomiyamaSwarm.C
    Making dependency list for source file noSwarm.C
    Making dependency list for source file newSwarmCorrection.C
    Making dependency list for source file swarmCorrection.C
    Making dependency list for source file IshiiZuber.C
    Making dependency list for source file WenYu.C
    Making dependency list for source file TomiyamaAnalytic.C
    Making dependency list for source file TomiyamaCorrelated.C
    Making dependency list for source file SyamlalOBrien.C
    Making dependency list for source file SchillerNaumann.C
    Making dependency list for source file Lain.C
    Making dependency list for source file GidaspowSchillerNaumann.C
    Making dependency list for source file GidaspowErgunWenYu.C
    Making dependency list for source file Gibilaro.C
    Making dependency list for source file Ergun.C
    Making dependency list for source file segregated.C
    Making dependency list for source file newDragModel.C
    Making dependency list for source file dragModel.C
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/dragModel/dragModel.C -o Make/linux64GccDPInt32Opt/dragModels/dragModel/dragModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/dragModel/newDragModel.C -o Make/linux64GccDPInt32Opt/dragModels/dragModel/newDragModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/segregated/segregated.C -o Make/linux64GccDPInt32Opt/dragModels/segregated/segregated.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/Ergun/Ergun.C -o Make/linux64GccDPInt32Opt/dragModels/Ergun/Ergun.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/Gibilaro/Gibilaro.C -o Make/linux64GccDPInt32Opt/dragModels/Gibilaro/Gibilaro.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C -o Make/linux64GccDPInt32Opt/dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C -o Make/linux64GccDPInt32Opt/dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/Lain/Lain.C -o Make/linux64GccDPInt32Opt/dragModels/Lain/Lain.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/SchillerNaumann/SchillerNaumann.C -o Make/linux64GccDPInt32Opt/dragModels/SchillerNaumann/SchillerNaumann.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/SyamlalOBrien/SyamlalOBrien.C -o Make/linux64GccDPInt32Opt/dragModels/SyamlalOBrien/SyamlalOBrien.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/TomiyamaCorrelated/TomiyamaCorrelated.C -o Make/linux64GccDPInt32Opt/dragModels/TomiyamaCorrelated/TomiyamaCorrelated.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/TomiyamaAnalytic/TomiyamaAnalytic.C -o Make/linux64GccDPInt32Opt/dragModels/TomiyamaAnalytic/TomiyamaAnalytic.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/WenYu/WenYu.C -o Make/linux64GccDPInt32Opt/dragModels/WenYu/WenYu.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c dragModels/IshiiZuber/IshiiZuber.C -o Make/linux64GccDPInt32Opt/dragModels/IshiiZuber/IshiiZuber.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c swarmCorrections/swarmCorrection/swarmCorrection.C -o Make/linux64GccDPInt32Opt/swarmCorrections/swarmCorrection/swarmCorrection.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c swarmCorrections/swarmCorrection/newSwarmCorrection.C -o Make/linux64GccDPInt32Opt/swarmCorrections/swarmCorrection/newSwarmCorrection.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c swarmCorrections/noSwarm/noSwarm.C -o Make/linux64GccDPInt32Opt/swarmCorrections/noSwarm/noSwarm.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c swarmCorrections/TomiyamaSwarm/TomiyamaSwarm.C -o Make/linux64GccDPInt32Opt/swarmCorrections/TomiyamaSwarm/TomiyamaSwarm.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c liftModels/liftModel/liftModel.C -o Make/linux64GccDPInt32Opt/liftModels/liftModel/liftModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c liftModels/liftModel/newLiftModel.C -o Make/linux64GccDPInt32Opt/liftModels/liftModel/newLiftModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c liftModels/noLift/noLift.C -o Make/linux64GccDPInt32Opt/liftModels/noLift/noLift.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c liftModels/constantLiftCoefficient/constantLiftCoefficient.C -o Make/linux64GccDPInt32Opt/liftModels/constantLiftCoefficient/constantLiftCoefficient.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c liftModels/Moraga/Moraga.C -o Make/linux64GccDPInt32Opt/liftModels/Moraga/Moraga.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c liftModels/LegendreMagnaudet/LegendreMagnaudet.C -o Make/linux64GccDPInt32Opt/liftModels/LegendreMagnaudet/LegendreMagnaudet.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c liftModels/TomiyamaLift/TomiyamaLift.C -o Make/linux64GccDPInt32Opt/liftModels/TomiyamaLift/TomiyamaLift.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c heatTransferModels/heatTransferModel/heatTransferModel.C -o Make/linux64GccDPInt32Opt/heatTransferModels/heatTransferModel/heatTransferModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c heatTransferModels/heatTransferModel/newHeatTransferModel.C -o Make/linux64GccDPInt32Opt/heatTransferModels/heatTransferModel/newHeatTransferModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c heatTransferModels/RanzMarshall/RanzMarshall.C -o Make/linux64GccDPInt32Opt/heatTransferModels/RanzMarshall/RanzMarshall.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c heatTransferModels/sphericalHeatTransfer/sphericalHeatTransfer.C -o Make/linux64GccDPInt32Opt/heatTransferModels/sphericalHeatTransfer/sphericalHeatTransfer.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c virtualMassModels/virtualMassModel/virtualMassModel.C -o Make/linux64GccDPInt32Opt/virtualMassModels/virtualMassModel/virtualMassModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c virtualMassModels/virtualMassModel/newVirtualMassModel.C -o Make/linux64GccDPInt32Opt/virtualMassModels/virtualMassModel/newVirtualMassModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c virtualMassModels/noVirtualMass/noVirtualMass.C -o Make/linux64GccDPInt32Opt/virtualMassModels/noVirtualMass/noVirtualMass.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c virtualMassModels/constantVirtualMassCoefficient/constantVirtualMassCoefficient.C -o Make/linux64GccDPInt32Opt/virtualMassModels/constantVirtualMassCoefficient/constantVirtualMassCoefficient.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c virtualMassModels/Lamb/Lamb.C -o Make/linux64GccDPInt32Opt/virtualMassModels/Lamb/Lamb.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c wallLubricationModels/wallLubricationModel/wallLubricationModel.C -o Make/linux64GccDPInt32Opt/wallLubricationModels/wallLubricationModel/wallLubricationModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c wallLubricationModels/wallLubricationModel/newWallLubricationModel.C -o Make/linux64GccDPInt32Opt/wallLubricationModels/wallLubricationModel/newWallLubricationModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c wallLubricationModels/noWallLubrication/noWallLubrication.C -o Make/linux64GccDPInt32Opt/wallLubricationModels/noWallLubrication/noWallLubrication.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c wallLubricationModels/Antal/Antal.C -o Make/linux64GccDPInt32Opt/wallLubricationModels/Antal/Antal.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c wallLubricationModels/Frank/Frank.C -o Make/linux64GccDPInt32Opt/wallLubricationModels/Frank/Frank.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c wallLubricationModels/TomiyamaWallLubrication/TomiyamaWallLubrication.C -o Make/linux64GccDPInt32Opt/wallLubricationModels/TomiyamaWallLubrication/TomiyamaWallLubrication.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.C -o Make/linux64GccDPInt32Opt/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c turbulentDispersionModels/turbulentDispersionModel/newTurbulentDispersionModel.C -o Make/linux64GccDPInt32Opt/turbulentDispersionModels/turbulentDispersionModel/newTurbulentDispersionModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c turbulentDispersionModels/noTurbulentDispersion/noTurbulentDispersion.C -o Make/linux64GccDPInt32Opt/turbulentDispersionModels/noTurbulentDispersion/noTurbulentDispersion.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c turbulentDispersionModels/constantTurbulentDispersionCoefficient/constantTurbulentDispersionCoefficient.C -o Make/linux64GccDPInt32Opt/turbulentDispersionModels/constantTurbulentDispersionCoefficient/constantTurbulentDispersionCoefficient.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c turbulentDispersionModels/Burns/Burns.C -o Make/linux64GccDPInt32Opt/turbulentDispersionModels/Burns/Burns.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c turbulentDispersionModels/Gosman/Gosman.C -o Make/linux64GccDPInt32Opt/turbulentDispersionModels/Gosman/Gosman.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c turbulentDispersionModels/LopezDeBertodano/LopezDeBertodano.C -o Make/linux64GccDPInt32Opt/turbulentDispersionModels/LopezDeBertodano/LopezDeBertodano.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c aspectRatioModels/aspectRatioModel/aspectRatioModel.C -o Make/linux64GccDPInt32Opt/aspectRatioModels/aspectRatioModel/aspectRatioModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c aspectRatioModels/aspectRatioModel/newAspectRatioModel.C -o Make/linux64GccDPInt32Opt/aspectRatioModels/aspectRatioModel/newAspectRatioModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c aspectRatioModels/constantAspectRatio/constantAspectRatio.C -o Make/linux64GccDPInt32Opt/aspectRatioModels/constantAspectRatio/constantAspectRatio.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c aspectRatioModels/TomiyamaAspectRatio/TomiyamaAspectRatio.C -o Make/linux64GccDPInt32Opt/aspectRatioModels/TomiyamaAspectRatio/TomiyamaAspectRatio.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c aspectRatioModels/VakhrushevEfremov/VakhrushevEfremov.C -o Make/linux64GccDPInt32Opt/aspectRatioModels/VakhrushevEfremov/VakhrushevEfremov.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c aspectRatioModels/Wellek/Wellek.C -o Make/linux64GccDPInt32Opt/aspectRatioModels/Wellek/Wellek.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c wallDependentModel/wallDependentModel.C -o Make/linux64GccDPInt32Opt/wallDependentModel/wallDependentModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -shared -Xlinker --add-needed -Xlinker --no-as-needed Make/linux64GccDPInt32Opt/dragModels/dragModel/dragModel.o Make/linux64GccDPInt32Opt/dragModels/dragModel/newDragModel.o Make/linux64GccDPInt32Opt/dragModels/segregated/segregated.o Make/linux64GccDPInt32Opt/dragModels/Ergun/Ergun.o Make/linux64GccDPInt32Opt/dragModels/Gibilaro/Gibilaro.o Make/linux64GccDPInt32Opt/dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.o Make/linux64GccDPInt32Opt/dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.o Make/linux64GccDPInt32Opt/dragModels/Lain/Lain.o Make/linux64GccDPInt32Opt/dragModels/SchillerNaumann/SchillerNaumann.o Make/linux64GccDPInt32Opt/dragModels/SyamlalOBrien/SyamlalOBrien.o Make/linux64GccDPInt32Opt/dragModels/TomiyamaCorrelated/TomiyamaCorrelated.o Make/linux64GccDPInt32Opt/dragModels/TomiyamaAnalytic/TomiyamaAnalytic.o Make/linux64GccDPInt32Opt/dragModels/WenYu/WenYu.o Make/linux64GccDPInt32Opt/dragModels/IshiiZuber/IshiiZuber.o Make/linux64GccDPInt32Opt/swarmCorrections/swarmCorrection/swarmCorrection.o Make/linux64GccDPInt32Opt/swarmCorrections/swarmCorrection/newSwarmCorrection.o Make/linux64GccDPInt32Opt/swarmCorrections/noSwarm/noSwarm.o Make/linux64GccDPInt32Opt/swarmCorrections/TomiyamaSwarm/TomiyamaSwarm.o Make/linux64GccDPInt32Opt/liftModels/liftModel/liftModel.o Make/linux64GccDPInt32Opt/liftModels/liftModel/newLiftModel.o Make/linux64GccDPInt32Opt/liftModels/noLift/noLift.o Make/linux64GccDPInt32Opt/liftModels/constantLiftCoefficient/constantLiftCoefficient.o Make/linux64GccDPInt32Opt/liftModels/Moraga/Moraga.o Make/linux64GccDPInt32Opt/liftModels/LegendreMagnaudet/LegendreMagnaudet.o Make/linux64GccDPInt32Opt/liftModels/TomiyamaLift/TomiyamaLift.o Make/linux64GccDPInt32Opt/heatTransferModels/heatTransferModel/heatTransferModel.o Make/linux64GccDPInt32Opt/heatTransferModels/heatTransferModel/newHeatTransferModel.o Make/linux64GccDPInt32Opt/heatTransferModels/RanzMarshall/RanzMarshall.o Make/linux64GccDPInt32Opt/heatTransferModels/sphericalHeatTransfer/sphericalHeatTransfer.o Make/linux64GccDPInt32Opt/virtualMassModels/virtualMassModel/virtualMassModel.o Make/linux64GccDPInt32Opt/virtualMassModels/virtualMassModel/newVirtualMassModel.o Make/linux64GccDPInt32Opt/virtualMassModels/noVirtualMass/noVirtualMass.o Make/linux64GccDPInt32Opt/virtualMassModels/constantVirtualMassCoefficient/constantVirtualMassCoefficient.o Make/linux64GccDPInt32Opt/virtualMassModels/Lamb/Lamb.o Make/linux64GccDPInt32Opt/wallLubricationModels/wallLubricationModel/wallLubricationModel.o Make/linux64GccDPInt32Opt/wallLubricationModels/wallLubricationModel/newWallLubricationModel.o Make/linux64GccDPInt32Opt/wallLubricationModels/noWallLubrication/noWallLubrication.o Make/linux64GccDPInt32Opt/wallLubricationModels/Antal/Antal.o Make/linux64GccDPInt32Opt/wallLubricationModels/Frank/Frank.o Make/linux64GccDPInt32Opt/wallLubricationModels/TomiyamaWallLubrication/TomiyamaWallLubrication.o Make/linux64GccDPInt32Opt/turbulentDispersionModels/turbulentDispersionModel/turbulentDispersionModel.o Make/linux64GccDPInt32Opt/turbulentDispersionModels/turbulentDispersionModel/newTurbulentDispersionModel.o Make/linux64GccDPInt32Opt/turbulentDispersionModels/noTurbulentDispersion/noTurbulentDispersion.o Make/linux64GccDPInt32Opt/turbulentDispersionModels/constantTurbulentDispersionCoefficient/constantTurbulentDispersionCoefficient.o Make/linux64GccDPInt32Opt/turbulentDispersionModels/Burns/Burns.o Make/linux64GccDPInt32Opt/turbulentDispersionModels/Gosman/Gosman.o Make/linux64GccDPInt32Opt/turbulentDispersionModels/LopezDeBertodano/LopezDeBertodano.o Make/linux64GccDPInt32Opt/aspectRatioModels/aspectRatioModel/aspectRatioModel.o Make/linux64GccDPInt32Opt/aspectRatioModels/aspectRatioModel/newAspectRatioModel.o Make/linux64GccDPInt32Opt/aspectRatioModels/constantAspectRatio/constantAspectRatio.o Make/linux64GccDPInt32Opt/aspectRatioModels/TomiyamaAspectRatio/TomiyamaAspectRatio.o Make/linux64GccDPInt32Opt/aspectRatioModels/VakhrushevEfremov/VakhrushevEfremov.o Make/linux64GccDPInt32Opt/aspectRatioModels/Wellek/Wellek.o Make/linux64GccDPInt32Opt/wallDependentModel/wallDependentModel.o -L/opt/openfoam5/platforms/linux64GccDPInt32Opt/lib \
        -L/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib -lcompressibleTwoPhaseSystemSrivastava -lcompressibleTransportModels -lfluidThermophysicalModels -lspecie  -o /home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libcompressibleEulerianInterfacialModelsSrivastava.so
    wmake phaseCompressibleTurbulenceModels
    wmakeLnInclude: linking include files to ./lnInclude
    Making dependency list for source file JohnsonJacksonParticleSlipFvPatchVectorField.C
    Making dependency list for source file JohnsonJacksonParticleThetaFvPatchScalarField.C
    Making dependency list for source file JohnsonJacksonSchaefferFrictionalStress.C
    Making dependency list for source file SchaefferFrictionalStress.C
    Making dependency list for source file JohnsonJacksonFrictionalStress.C
    Making dependency list for source file newFrictionalStressModel.C
    Making dependency list for source file frictionalStressModel.C
    Making dependency list for source file SyamlalRogersOBrienPressure.C
    Making dependency list for source file LunPressure.C
    Making dependency list for source file newGranularPressureModel.C
    Making dependency list for source file granularPressureModel.C
    Making dependency list for source file SinclairJacksonRadial.C
    Making dependency list for source file LunSavageRadial.C
    Making dependency list for source file CarnahanStarlingRadial.C
    Making dependency list for source file newRadialModel.C
    Making dependency list for source file radialModel.C
    Making dependency list for source file HrenyaSinclairConductivity.C
    Making dependency list for source file SyamlalConductivity.C
    Making dependency list for source file GidaspowConductivity.C
    Making dependency list for source file newConductivityModel.C
    Making dependency list for source file conductivityModel.C
    Making dependency list for source file noneViscosity.C
    Making dependency list for source file HrenyaSinclairViscosity.C
    Making dependency list for source file SyamlalViscosity.C
    Making dependency list for source file GidaspowViscosity.C
    Making dependency list for source file newViscosityModel.C
    Making dependency list for source file viscosityModel.C
    Making dependency list for source file kineticTheoryModel.C
    Making dependency list for source file phasePressureModel.C
    Making dependency list for source file phaseCompressibleTurbulenceModels.C
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c phaseCompressibleTurbulenceModels.C -o Make/linux64GccDPInt32Opt/phaseCompressibleTurbulenceModels.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c phasePressureModel/phasePressureModel.C -o Make/linux64GccDPInt32Opt/phasePressureModel/phasePressureModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/kineticTheoryModel/kineticTheoryModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/kineticTheoryModel/kineticTheoryModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/viscosityModel/viscosityModel/viscosityModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/viscosityModel/viscosityModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/viscosityModel/viscosityModel/newViscosityModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/viscosityModel/newViscosityModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/viscosityModel/Gidaspow/GidaspowViscosity.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/Gidaspow/GidaspowViscosity.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/viscosityModel/Syamlal/SyamlalViscosity.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/Syamlal/SyamlalViscosity.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/viscosityModel/none/noneViscosity.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/none/noneViscosity.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/conductivityModel/conductivityModel/conductivityModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/conductivityModel/conductivityModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/conductivityModel/conductivityModel/newConductivityModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/conductivityModel/newConductivityModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/conductivityModel/Gidaspow/GidaspowConductivity.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/Gidaspow/GidaspowConductivity.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/conductivityModel/Syamlal/SyamlalConductivity.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/Syamlal/SyamlalConductivity.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/radialModel/radialModel/radialModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/radialModel/radialModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/radialModel/radialModel/newRadialModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/radialModel/newRadialModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/radialModel/CarnahanStarling/CarnahanStarlingRadial.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/CarnahanStarling/CarnahanStarlingRadial.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/radialModel/LunSavage/LunSavageRadial.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/LunSavage/LunSavageRadial.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/radialModel/SinclairJackson/SinclairJacksonRadial.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/SinclairJackson/SinclairJacksonRadial.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/granularPressureModel/granularPressureModel/granularPressureModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/granularPressureModel/granularPressureModel/granularPressureModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/granularPressureModel/granularPressureModel/newGranularPressureModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/granularPressureModel/granularPressureModel/newGranularPressureModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/granularPressureModel/Lun/LunPressure.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/granularPressureModel/Lun/LunPressure.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/frictionalStressModel/frictionalStressModel/frictionalStressModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/frictionalStressModel/frictionalStressModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/frictionalStressModel/newFrictionalStressModel.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/Schaeffer/SchaefferFrictionalStress.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleTheta/JohnsonJacksonParticleThetaFvPatchScalarField.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleTheta/JohnsonJacksonParticleThetaFvPatchScalarField.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.C -o Make/linux64GccDPInt32Opt/kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/transportModel -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -shared -Xlinker --add-needed -Xlinker --no-as-needed Make/linux64GccDPInt32Opt/phaseCompressibleTurbulenceModels.o Make/linux64GccDPInt32Opt/phasePressureModel/phasePressureModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/kineticTheoryModel/kineticTheoryModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/viscosityModel/viscosityModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/viscosityModel/newViscosityModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/Gidaspow/GidaspowViscosity.o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/Syamlal/SyamlalViscosity.o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.o Make/linux64GccDPInt32Opt/kineticTheoryModels/viscosityModel/none/noneViscosity.o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/conductivityModel/conductivityModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/conductivityModel/newConductivityModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/Gidaspow/GidaspowConductivity.o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/Syamlal/SyamlalConductivity.o Make/linux64GccDPInt32Opt/kineticTheoryModels/conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/radialModel/radialModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/radialModel/newRadialModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/CarnahanStarling/CarnahanStarlingRadial.o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/LunSavage/LunSavageRadial.o Make/linux64GccDPInt32Opt/kineticTheoryModels/radialModel/SinclairJackson/SinclairJacksonRadial.o Make/linux64GccDPInt32Opt/kineticTheoryModels/granularPressureModel/granularPressureModel/granularPressureModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/granularPressureModel/granularPressureModel/newGranularPressureModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/granularPressureModel/Lun/LunPressure.o Make/linux64GccDPInt32Opt/kineticTheoryModels/granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/frictionalStressModel/frictionalStressModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/frictionalStressModel/newFrictionalStressModel.o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/Schaeffer/SchaefferFrictionalStress.o Make/linux64GccDPInt32Opt/kineticTheoryModels/frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.o Make/linux64GccDPInt32Opt/kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleTheta/JohnsonJacksonParticleThetaFvPatchScalarField.o Make/linux64GccDPInt32Opt/kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.o -L/opt/openfoam5/platforms/linux64GccDPInt32Opt/lib \
        -L/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib -lcompressibleTransportModels -lfluidThermophysicalModels -lspecie -lturbulenceModels -lcompressibleTurbulenceModels -lincompressibleTransportModels -lcompressibleTwoPhaseSystemSrivastava -lcompressibleEulerianInterfacialModelsSrivastava -lfiniteVolume -lfvOptions -lmeshTools  -o /home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsSrivastava.so
    Making dependency list for source file twoPhaseEulerSrivastavaFoam.C
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -c twoPhaseEulerSrivastavaFoam.C -o Make/linux64GccDPInt32Opt/twoPhaseEulerSrivastavaFoam.o
    g++ -std=c++11 -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -O3  -DNoRepository -ftemplate-depth-100 -I/opt/openfoam5/src/transportModels/compressible/lnInclude -I/opt/openfoam5/src/thermophysicalModels/basic/lnInclude -I/opt/openfoam5/src/transportModels/incompressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/turbulenceModels/lnInclude -I/opt/openfoam5/src/TurbulenceModels/compressible/lnInclude -I/opt/openfoam5/src/TurbulenceModels/phaseCompressible/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/interfacialModels/lnInclude -I/home/apr207/OpenFOAM/apr207-5.0/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/twoPhaseSystem/lnInclude -I/opt/openfoam5/src/finiteVolume/lnInclude -I/opt/openfoam5/src/meshTools/lnInclude -I/opt/openfoam5/src/sampling/lnInclude -IlnInclude -I. -I/opt/openfoam5/src/OpenFOAM/lnInclude -I/opt/openfoam5/src/OSspecific/POSIX/lnInclude   -fPIC -Xlinker --add-needed -Xlinker --no-as-needed Make/linux64GccDPInt32Opt/twoPhaseEulerSrivastavaFoam.o -L/opt/openfoam5/platforms/linux64GccDPInt32Opt/lib \
        -L/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib -lcompressibleTransportModels -lfluidThermophysicalModels -lspecie -lturbulenceModels -lcompressibleTurbulenceModels -lphaseCompressibleTurbulenceModelsSrivastava -lincompressibleTransportModels -lcompressibleTwoPhaseSystemSrivastava -lcompressibleEulerianInterfacialModelsSrivastava -lfiniteVolume -lfvOptions -lmeshTools -lsampling -lOpenFOAM -ldl  \
        -lm -o /home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/bin/twoPhaseEulerSrivastavaFoam

        

        
Add new friction model ReugeRahimi
----------------------------------

Simply copy existing JohnsonJacksonSchaeffer and see if it works:

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/kineticTheoryModels/frictionalStressModel
    mkdir ReugeRahimi
    cd ReugeRahimi
    cd ..
    cp -rf JohnsonJacksonSchaeffer/*.*  ReugeRahimi/.
    
    cd ReugeRahimi

    mv JohnsonJacksonSchaefferFrictionalStress.C ReugeRahimiFrictionalStress.C
    mv JohnsonJacksonSchaefferFrictionalStress.H ReugeRahimiFrictionalStress.H
    sed -i s/JohnsonJacksonSchaeffer/ReugeRahimi/g ReugeRahimiFrictionalStress.C
    sed -i s/JohnsonJacksonSchaeffer/ReugeRahimi/g ReugeRahimiFrictionalStress.H  

Added the line:    
    
::
    
    cd ../../..
    nano Make/files
        kineticTheoryModels/frictionalStressModel/ReugeRahimi/ReugeRahimiFrictionalStress.C
        
Recompile:        
        
::

    ./Allwclean
    ./Allwmake

    
Now add the arguments Theta and da for the new friction model to the nu function:

::

    frictionalStressModel.H
    JohnsonJacksonFrictionalStress.H and .C
    SchaefferFrictionalStress.H and .C
    ReugeRahimiFrictionalStress.H and .C
    JohnsonJacksonFrictionalStress.H and .C
    kineticTheoryModel.C

Add these, e.g. for the JohnsonJackson case:    
    
::

    Foam::tmp<Foam::volScalarField>
    Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::nu
    (
        const phaseModel& phase,
        const dimensionedScalar& alphaMinFriction,
        const dimensionedScalar& alphaMax,
        const volScalarField& pf,
        const volSymmTensorField& D,
        const volScalarField& Theta,
        const dimensionedScalar& da
    ) const
    {
        return dimensionedScalar("0.5", dimTime, 0.5)*pf*sin(phi_);
    }  

Also add Theta and da to the call in kineticTheoryModel.C    

::

    nuFric_ = frictionalStressModel_->nu
    (
        phase_,
        alphaMinFriction_,
        alphaMax_,
        pf/rho,
        D,
        Theta_,
        da
    );
    
    
    
Now check recompilation

::

    ./Allwclean
    ./Allwmake

Boundary correction looks to be inconsitent with rest of friction model - Add new friction model ReugeRahimiBoundary
--------------------------------------------------------------------------------------------------------------------

Simply copy existing JohnsonJacksonSchaeffer and see if it works:

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/kineticTheoryModels/frictionalStressModel
    mkdir ReugeRahimiBoundary
    cd ReugeRahimiBoundary
    cd ..
    cp -rf ReugeRahimi/*.*  ReugeRahimiBoundary/.
    
    cd ReugeRahimiBoundary

    mv ReugeRahimiFrictionalStress.C ReugeRahimiBoundaryFrictionalStress.C
    mv ReugeRahimiFrictionalStress.H ReugeRahimiBoundaryFrictionalStress.H
    sed -i s/ReugeRahimi/ReugeRahimiBoundary/g ReugeRahimiBoundaryFrictionalStress.C
    sed -i s/ReugeRahimi/ReugeRahimiBoundary/g ReugeRahimiBoundaryFrictionalStress.H  

Added the line:    
    
::
    
    cd ../../..
    nano Make/files
        kineticTheoryModels/frictionalStressModel/ReugeRahimiBoundary/ReugeRahimiBoundaryFrictionalStress.C
        
Recompile:        
        
::

    ./Allwclean
    ./Allwmake

    

 Copy to emps-kabalevsky
 -----------------------
 
::
 
    five
    cd $WM_PROJECT_USER_DIR
    scp -r applications apr207@emps-kabalevsky:/home/links/apr207/OpenFOAM/apr207-5.0  
    
::

    ./Allwclean
    ./Allwmake
    
Now just modify the friction model and recompile using:

::

    wclean libso phaseCompressibleTurbulenceModels
    wmake phaseCompressibleTurbulenceModels
        
        
Add Srivastava Friction Model without critical state hypothesis 
---------------------------------------------------------------        

Simply copy existing ReugeRahimiBoundary and see if it works:

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/kineticTheoryModels/frictionalStressModel
    mkdir Srivastava
    cp -rf ReugeRahimiBoundary/*.*  Srivastava/.
    
    cd Srivastava

    mv ReugeRahimiBoundaryFrictionalStress.C SrivastavaFrictionalStress.C
    mv ReugeRahimiBoundaryFrictionalStress.H SrivastavaFrictionalStress.H
    sed -i s/ReugeRahimiBoundary/Srivastava/g SrivastavaFrictionalStress.C
    sed -i s/ReugeRahimiBoundary/Srivastava/g SrivastavaFrictionalStress.H  

Added the line:    
    
::
    
    cd ../../..
    nano Make/files
        kineticTheoryModels/frictionalStressModel/Srivastava/SrivastavaFrictionalStress.C
        
Recompile:        
        
::

    ./Allwclean
    ./Allwmake


        
Add ReugeRahimiCoulomb
----------------------

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/kineticTheoryModels/frictionalStressModel
    mkdir ReugeRahimiCoulomb
    cp -rf ReugeRahimi/*.*  ReugeRahimiCoulomb/.
    
    cd ReugeRahimiCoulomb

    mv ReugeRahimiFrictionalStress.C ReugeRahimiCoulombFrictionalStress.C
    mv ReugeRahimiFrictionalStress.H ReugeRahimiCoulombFrictionalStress.H
    sed -i s/ReugeRahimi/ReugeRahimiCoulomb/g ReugeRahimiCoulombFrictionalStress.C
    sed -i s/ReugeRahimi/ReugeRahimiCoulomb/g ReugeRahimiCoulombFrictionalStress.H  

Added the line:    
    
::
    
    cd ../../..
    nano Make/files
        kineticTheoryModels/frictionalStressModel/ReugeRahimiCoulomb/ReugeRahimiCoulombFrictionalStress.C
        
Recompile:        
        
::

    ./Allwclean
    ./Allwmake

Add new viscosity model - Srivastava
------------------------------------

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/kineticTheoryModels/viscosityModel
    mkdir Srivastava
    cp -rf Syamlal/*.*  Srivastava/.
    
    cd Srivastava

    mv SyamlalViscosity.C SrivastavaViscosity.C
    mv SyamlalViscosity.H SrivastavaViscosity.H
    sed -i s/Syamlal/Srivastava/g SrivastavaViscosity.C
    sed -i s/Syamlal/Srivastava/g SrivastavaViscosity.H  

Added the line:    
    
::
    
    cd ../../..
    nano Make/files
        kineticTheoryModels/viscosityModel/Srivastava/SrivastavaViscosity.C
        
Recompile:        
        
::
    cd ..
    ./Allwclean
    ./Allwmake

Add new conductivity model - Srivastava
---------------------------------------

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/kineticTheoryModels/conductivityModel
    mkdir Srivastava
    cp -rf Syamlal/*.*  Srivastava/.
    
    cd Srivastava

    mv SyamlalConductivity.C SrivastavaConductivity.C
    mv SyamlalConductivity.H SrivastavaConductivity.H
    sed -i s/Syamlal/Srivastava/g SrivastavaConductivity.C
    sed -i s/Syamlal/Srivastava/g SrivastavaConductivity.H  

Added the line:    
    
::
    
    cd ../../..
    nano Make/files
        kineticTheoryModels/conductivityModel/Srivastava/SrivastavaConductivity.C
        
Recompile:        
        
::
    cd ..
    ./Allwclean
    ./Allwmake    
    
Add new Friction Model - SrivastavaCriticalStateHypothesis
----------------------------------------------------------

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/twoPhaseEulerSrivastavaFoam/phaseCompressibleTurbulenceModels/kineticTheoryModels/frictionalStressModel
    mkdir SrivastavaCriticalStateHypothesis
    cp -rf ReugeRahimiCoulomb/*.*  SrivastavaCriticalStateHypothesis/.
    
    cd SrivastavaCriticalStateHypothesis

    mv ReugeRahimiCoulombFrictionalStress.C SrivastavaCriticalStateHypothesisFrictionalStress.C
    mv ReugeRahimiCoulombFrictionalStress.H SrivastavaCriticalStateHypothesisFrictionalStress.H
    sed -i s/ReugeRahimiCoulomb/SrivastavaCriticalStateHypothesis/g SrivastavaCriticalStateHypothesisFrictionalStress.C
    sed -i s/ReugeRahimiCoulomb/SrivastavaCriticalStateHypothesis/g SrivastavaCriticalStateHypothesisFrictionalStress.H  

Added the line:    
    
::
    
    cd ../../..
    nano Make/files
        kineticTheoryModels/frictionalStressModel/SrivastavaCriticalStateHypothesis/SrivastavaCriticalStateHypothesisFrictionalStress.C
        
Recompile:        
        
::

    ./Allwclean
    ./Allwmake
    
    
    
    
    
    
    
    
    
    
    
Create new solver - twoPhaseEulerFrictionFoam
---------------------------------------------
 
 
 Copy the code from source

::

    five
    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase
    cp -rf twoPhaseEulerSrivastavaFoam twoPhaseEulerFrictionFoam
    cd twoPhaseEulerFrictionFoam
    mv twoPhaseEulerSrivastavaFoam.C twoPhaseEulerFrictionFoam.C
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerFrictionFoam/g twoPhaseEulerFrictionFoam.C
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerFrictionFoam/g Make/files
    ./Allwclean
    ./Allwmake

    
    
Correct Make and options files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

twoPhaseEulerFrictionFoam    
""""""""""""""""""""""""" 

nano Make/files (altered) - location and name of executable has changed

The first line means the .C file to compile is located one directory above the Make folder.

The second line means put the executable in the directory $FOAM_USER_APPBIN

::

    twoPhaseEulerFrictionFoam.C

    EXE = $(FOAM_USER_APPBIN)/twoPhaseEulerFrictionFoam

https://cfd.direct/openfoam/user-guide/compiling-applications/    
    
nano Make/options (altered):

* EXE_INC is the path to the libraries
* EXE_LIBS is the library names (these have been altered). The library path must also be given, the user location.
* The actual library files to be linked must be specified using the -l option and removing the lib prefix and .so extension from the library file name, e.g.  libnew.so is included with the flag -lnew

Why is -L specified only once? Possibly because -L adds an additional path to the default path ($FOAM_LIBBIN)

::

    EXE_INC = \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/phaseCompressibleTurbulenceModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/interfacialModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/twoPhaseSystem/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/sampling/lnInclude

    EXE_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie \
        -lturbulenceModels \
        -lcompressibleTurbulenceModels \
        -lphaseCompressibleTurbulenceModelsFriction \
        -lincompressibleTransportModels \
        -lcompressibleTwoPhaseSystemFriction \
        -lcompressibleEulerianInterfacialModelsFriction \
        -lfiniteVolume \
        -lfvOptions \
        -lmeshTools \
        -lsampling

The result:



        
        
interfacialModels    
"""""""""""""""""             

The location of the nine .C files to compile is the same relative to the Make folder.

nano interfacialModels/Make/files (altered):

* The location and name of the EulerianInterfacialModels library has changed

::

    LIB = $(FOAM_USER_LIBBIN)/libcompressibleEulerianInterfacialModelsFriction

nano interfacialModels/Make/options (altered):

* The name and location of the phaseModel library has changed
    
::

    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/transportModel \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/twoPhaseSystem/lnInclude

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lcompressibleTwoPhaseSystemFriction \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie

        
phaseCompressibleTurbulenceModels
"""""""""""""""""""""""""""""""""

nano phaseCompressibleTurbulenceModels/Make/files (altered):

::

    LIB = $(FOAM_USER_LIBBIN)/libphaseCompressibleTurbulenceModelsFriction
       
        
nano phaseCompressibleTurbulenceModels/Make/options (altered):        

::

    EXE_INC = \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/twoPhaseSystem/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/interfacialModels/lnInclude\
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/transportModel \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie \
        -lturbulenceModels \
        -lcompressibleTurbulenceModels \
        -lincompressibleTransportModels \
        -lcompressibleTwoPhaseSystemFriction \
        -lcompressibleEulerianInterfacialModelsFriction \
        -lfiniteVolume \
        -lfvOptions \
        -lmeshTools
        
        
twoPhaseSystem
""""""""""""""

nano twoPhaseSystem/Make/files

::

     LIB = $(FOAM_USER_LIBBIN)/libcompressibleTwoPhaseSystemFriction   
  
  
nano twoPhaseSystem/Make/options  
  
::

    EXE_INC = \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/twoPhaseSystem \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/interfacialModels/lnInclude \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/sampling/lnInclude

    LIB_LIBS = \
        -lincompressibleTransportModels \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie

Recompile:        
        
::

    ./Allwclean
    ./Allwmake
    
Remove friction from kineticTheory.C
------------------------------------

Remove frictionalPressurePrime from pPrime

::

    Foam::tmp<Foam::volScalarField>
    Foam::RASModels::kineticTheoryModel::pPrime() const
    {
        const volScalarField& rho = phase_.rho();

        tmp<volScalarField> tpPrime
        (
            Theta_
        *granularPressureModel_->granularPressureCoeffPrime
            (
                alpha_,
                radialModel_->g0(alpha_, alphaMinFriction_, alphaMax_),
                radialModel_->g0prime(alpha_, alphaMinFriction_, alphaMax_),
                rho,
                e_
            )
        /*
        +  frictionalStressModel_->frictionalPressurePrime
            (
                phase_,
                alphaMinFriction_,
                alphaMax_
            )
        */
        );

        volScalarField::Boundary& bpPrime =
            tpPrime.ref().boundaryFieldRef();

        forAll(bpPrime, patchi)
        {
            if (!bpPrime[patchi].coupled())
            {
                bpPrime[patchi] == 0;
            }
        }

        return tpPrime;
    }


Remove nuFric from nut_

::

    //nut_ += nuFric_; 



Create new RASModel called frictionModel
----------------------------------------

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/phaseCompressibleTurbulenceModels
    mkdir frictionModel
    cp -rf kineticTheory/ frictionModel
    cd frictionModel
    mv kineticTheoryModel frictionModel
    rm -rf conductivityModel derivedFvPatchFields granularPressureModel radialModel viscosityModel
    cd frictionModel
    mv kineticTheoryModel.C frictionModel.C
    mv kineticTheoryModel.H frictionModel.H
    sed -i s/kineticTheoryModel/frictionModel/g frictionModel.C
    sed -i s/kineticTheoryModel/frictionModel/g frictionModel.H
    sed -i s/kineticTheoryModel/frictionModel/g Make/files

Ad hoc modiciations    
    
::

    TypeName("frictionModel");
    replace fricionModels with frictionModel
    
    
replace kineticTheoryModels with frictionModel

::

    cd frictionalStressModel/frictionalStressModel
    sed -i s/kineticTheoryModels/frictionModel/g frictionalStressModel.C
    sed -i s/kineticTheoryModels/frictionModel/g frictionalStressModel.H
    sed -i s/kineticTheoryModels/frictionModel/g newFrictionalStressModel.C
    
::

    cd ../SrivastavaCriticalStateHypothesis
    sed -i s/kineticTheoryModels/frictionModel/g SrivastavaCriticalStateHypothesisFrictionalStress.C
    sed -i s/kineticTheoryModels/frictionModel/g SrivastavaCriticalStateHypothesisFrictionalStress.H
    
Add these lines to phaseCompressibleTurbulenceModels/Make/files
    
::

    frictionModel/frictionModel/frictionModel.C
    frictionModel/frictionalStressModel/SrivastavaCriticalStateHypothesis/SrivastavaCriticalStateHypothesisFrictionalStress.C
    
    
    
::

    sed -i s/SrivastavaCriticalStateHypothesis/SrivastavaSundaresan/g SrivastavaSundaresan.C
    sed -i s/SrivastavaCriticalStateHypothesis/SrivastavaSundaresan/g SrivastavaSundaresan.H    
    

    
    
    
    
::

    sed -i s/frictionalStressModel/frictionalStress/g frictionalStress.C
    sed -i s/frictionalStressModel/frictionalStress/g frictionalStress.H    
    sed -i s/frictionalStressModel/frictionalStress/g newFrictionalStress.C   
    
    
::

    frictionModel/frictionModel/frictionModel.C
    frictionModel/frictionalStress/SrivastavaSundaresan/SrivastavaSundaresan.C

::

    sed -i s/frictionalStressModel/frictionalStress/g frictionModel.C
    sed -i s/frictionalStressModel/frictionalStress/g frictionModel.H    
    
    
    
    
Create New Boundary Condition
-----------------------------

::

    cp -rf JohnsonJacksonParticleTheta JohnsonJacksonParticleVelocity

    cd JohnsonJacksonParticleVelocity

    mv JohnsonJacksonParticleThetaFvPatchScalarField.C JohnsonJacksonParticleVelocityFvPatchScalarField.C
    mv JohnsonJacksonParticleThetaFvPatchScalarField.H JohnsonJacksonParticleVelocityFvPatchScalarField.H

    sed -i s/JohnsonJacksonParticleTheta/JohnsonJacksonParticleVelocity/g JohnsonJacksonParticleVelocityFvPatchScalarField.C
    sed -i s/JohnsonJacksonParticleTheta/JohnsonJacksonParticleVelocity/g JohnsonJacksonParticleVelocityFvPatchScalarField.H


phaseCompressibleTurbulenceModels/Make/files

::

    kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleVelocity/JohnsonJacksonParticleVelocityFvPatchScalarField.C


Need to change FvPatchScalarField to FvPatchVectorField


::

    mv JohnsonJacksonParticleVelocityFvPatchScalarField.C JohnsonJacksonParticleVelocityFvPatchVectorField.C
    mv JohnsonJacksonParticleVelocityFvPatchScalarField.H JohnsonJacksonParticleVelocityFvPatchVectorField.H
    
    sed -i s/ScalarField/ScalarVector/g JohnsonJacksonParticleVelocityFvPatchVectorField.C
    sed -i s/ScalarField/ScalarVector/g JohnsonJacksonParticleVelocityFvPatchVectorField.H
    
    sed -i s/ScalarVector/VectorField/g JohnsonJacksonParticleVelocityFvPatchVectorField.C
    sed -i s/ScalarVector/VectorField/g JohnsonJacksonParticleVelocityFvPatchVectorField.H 
    
::

    kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleVelocity/JohnsonJacksonParticleVelocityFvPatchVectorField.C    
    
Replace scalar with vector

::
    
    sed -i s/scalar/vector/g JohnsonJacksonParticleVelocityFvPatchVectorField.C
    sed -i s/scalar/vector/g JohnsonJacksonParticleVelocityFvPatchVectorField.H 
    
Send errors to log file!   

::

    ./Allwmake > logmake 2>&1


Create New Boundary Condition based on JohnsonJacksonParticleSlip
-----------------------------------------------------------------

::

    cp -rf JohnsonJacksonParticleSlip JohnsonJacksonParticleSlipMixed

    cd JohnsonJacksonParticleSlipMixed

    mv JohnsonJacksonParticleSlipFvPatchVectorField.C JohnsonJacksonParticleSlipMixedFvPatchVectorField.C
    mv JohnsonJacksonParticleSlipFvPatchVectorField.H JohnsonJacksonParticleSlipMixedFvPatchVectorField.H

    sed -i s/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipMixed/g JohnsonJacksonParticleSlipMixedFvPatchVectorField.C
    sed -i s/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipMixed/g JohnsonJacksonParticleSlipMixedFvPatchVectorField.H


phaseCompressibleTurbulenceModels/Make/files

::

    kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleSlipMixed/JohnsonJacksonParticleSlipMixedFvPatchVectorField.C

    
Send errors to log file! - This works   

::

    ./Allwmake > logmake 2>&1    
    
    
    
need to replace partialSlip with mixed 

::

    sed -i s/partialSlipFvPatchVectorField/mixedFvPatchVectorField/g JohnsonJacksonParticleSlipMixedFvPatchVectorField.C
    sed -i s/partialSlipFvPatchVectorField/mixedFvPatchVectorField/g JohnsonJacksonParticleSlipMixedFvPatchVectorField.H


    
in JohnsonJacksonParticleSlipMixedFvPatchVectorField.H
    
::

    #include "partialSlipFvPatchFields.H" #include "mixedFvPatchFields.H" 
    
    
now add angleOfFriction



Copy twoPhaseEulerSrivastavaFoam to twoPhaseEulerSrivastavaFrictionFoam
-----------------------------------------------------------------------

files

::

    cp -rf twoPhaseEulerSrivastavaFoam twoPhaseEulerSrivastavaFrictionFoam
    cd twoPhaseEulerSrivastavaFrictionFoam
    mv twoPhaseEulerSrivastavaFoam.C twoPhaseEulerSrivastavaFrictionFoam.C
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g twoPhaseEulerSrivastavaFrictionFoam.C
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g Make/options
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g phaseCompressibleTurbulenceModels/Make/options
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g interfacialModels/Make/options    
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaFrictionFoam/g twoPhaseSystem/Make/options
    ./Allwclean
    ./Allwmake
    

copy frictionModel 

::

    cp -rf twoPhaseEulerFrictionFoam/frictionModels/ twoPhaseEulerSrivastavaFrictionFoam/frictionModels
    sed -i s/twoPhaseEulerFrictionFoam/twoPhaseEulerSrivastavaFrictionFoam/g frictionModels/Make/files
    sed -i s/twoPhaseEulerFrictionFoam/twoPhaseEulerSrivastavaFrictionFoam/g frictionModels/Make/options
    
    
add friction model to Make/options

::

    -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFrictionFoam/frictionModels/lnInclude \
    

Remove friction:

kineticTheoryModel.C

::
    
    Foam::tmp<Foam::volScalarField>
    Foam::RASModels::kineticTheoryModel::pPrime() const
    {
        const volScalarField& rho = phase_.rho();

        tmp<volScalarField> tpPrime
        (
            Theta_
        *granularPressureModel_->granularPressureCoeffPrime
            (
                alpha_,
                radialModel_->g0(alpha_, alphaMinFriction_, alphaMax_),
                radialModel_->g0prime(alpha_, alphaMinFriction_, alphaMax_),
                rho,
                e_
            )
    //     +  frictionalStressModel_->frictionalPressurePrime
    //        (
    //            phase_,
    //            alphaMinFriction_,
    //            alphaMax_
    //        )
        );

        volScalarField::Boundary& bpPrime =
            tpPrime.ref().boundaryFieldRef();

        forAll(bpPrime, patchi)
        {
            if (!bpPrime[patchi].coupled())
            {
                bpPrime[patchi] == 0;
            }
        }

        return tpPrime;
    }    
    
    
::

    // pf_ =  frictionalStressModel_->frictionalPressure
    // (
    //    phase_,
    //    alphaMinFriction_,
    //    alphaMax_
    // );   
    
    
    //    Info<< "max(pf) = " << max(pf_).value() << endl;
    //    Info<< "min(pf) = " << min(pf_).value() << endl;

    //    nuFric_ = frictionalStressModel_->nu
    //    (
    //        phase_,
    //        alphaMinFriction_,
    //        alphaMax_,
    //        pf_/rho,
    //        Up,
    //        divUp,
    //        D,
    //        Theta_,
    //        da,
    //        wallFriction_
    //    );

        // Limit viscosity and add frictional viscosity
        nut_.min(maxNut_);
        // nuFric_ = min(nuFric_, maxNut_ - nut_);
        nuKinColl_ = nut_;
        // nut_ += nuFric_;
    }
    
    Info<< "max(nut) = " << max(nut_).value() << endl;
    Info<< "max(nuKinColl) = " << max(nuKinColl_).value() << endl;
    // Info<< "max(nuFric) = " << max(nuFric_).value() << endl;
    
    
Need to associate the new friction model with the phaseModel object e.g.

::

    phase1.sigmaFIso(phase1, alphaMinFriction, alphaMax)

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    volScalarField Theta = U.mesh().lookupObject<volScalarField>("Theta");

    phase1.sigmaFDev(phase1, alphaMinFriction, alphaMax, sigmaFIso, D, Theta, da, wallFriction)

Create frictionSourceModel as sub-model of twoPhaseSystem.

::

    cp -rf diameterModels frictionSourceModels
    cd frictionSourceModels
    mv diameterModel frictionSourceModel
    cd frictionSourceModel
    mv diameterModel.C frictionSourceModel.C
    mv diameterModel.H frictionSourceModel.H
    mv newDiameterModel.C newFrictionSourceModel.C
    
Replace words:

::

    sed -i s/diameterModel/frictionSourceModel/g frictionSourceModel.H
    sed -i s/diameterModel/frictionSourceModel/g frictionSourceModel.C
    sed -i s/diameterModel/frictionSourceModel/g newFrictionSourceModel.C
    

Remove everything except constantDiameter:

::

    rm -rf IATE
    rm -rf isothermalDiameter
    

Rename 

::

    mv constantDiameter.H constantFriction.H
    mv constantDiameter.C constantFriction.C
    sed -i s/diameterModel/frictionSourceModel/g constantFriction.H
    sed -i s/diameterModel/frictionSourceModel/g constantFriction.C
    sed -i s/constantDiameter/constantFriction/g constantFriction.C
    
Add to Make files:

::

    frictionSourceModels/frictionSourceModel/frictionSourceModel.C
    frictionSourceModels/frictionSourceModel/newFrictionSourceModel.C
    frictionSourceModels/constantFriction/constantFriction.C


Replace function d() with function frictionalStressIso()    
    
    
    
Got to here:

    
    
    
New friction model:    
    
::

    frictionPtr_ = frictionSourceModel::New
    (
        phaseDict_,
        *this
    );



Change constantFriction to SrivastavaFriction

::

    mv constantFriction.C SrivastavaFriction.C
    mv constantFriction.H SrivastavaFriction.H
    
::

    sed -i s/constant/SrivastavaFriction/g SrivastavaFriction.H
    sed -i s/constant/SrivastavaFriction/g SrivastavaFriction.C
    
    
Change frictionSourceModel to frictionModel

::

    sed -i s/frictionSourceModel/frictionModel/g frictionModel.H
    sed -i s/frictionSourceModel/frictionModel/g frictionModel.C

Change SrivastavaFriction to Srivastava

::

    sed -i s/SrivastavaFriction/Srivastava/g Srivastava.H
    sed -i s/SrivastavaFriction/Srivastava/g Srivastava.C
    sed -i s/frictionSourceModel/frictionModel/g Srivastava.H
    sed -i s/frictionSourceModel/frictionModel/g Srivastava.C

Remove old friction model - Make/files

::

    -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerSrivastavaFrictionFoam/frictionModels/lnInclude \

    
diameterProperties becomes frictionProperties    
    
::

    sed -i s/diameterProperties/frictionProperties/g frictionModel.H
    sed -i s/diameterProperties/frictionProperties/g frictionModel.C    
    sed -i s/diameterProperties/frictionProperties/g newFrictionModel.C
   
   
   
::

    sed -i s/diameterProperties/frictionProperties/g Srivastava.C
    sed -i s/diameterProperties/frictionProperties/g Srivastava.H     
    
Remove links to frictionalStressModel
    
::
    
    kineticTheoryModels/frictionalStressModel/frictionalStressModel/frictionalStressModel.C
    kineticTheoryModels/frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C
    kineticTheoryModels/frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/ReugeRahimi/ReugeRahimiFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/ReugeRahimiBoundary/ReugeRahimiBoundaryFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/Srivastava/SrivastavaFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/ReugeRahimiCoulomb/ReugeRahimiCoulombFrictionalStress.C
    kineticTheoryModels/frictionalStressModel/SrivastavaCriticalStateHypothesis/SrivastavaCriticalStateHypothesisFrictionalStress.C
    
    
remove frictionalStress.H from kineticTheoryModel.H   
also remove all references to frictionalStress object


Need to add None for the frictionModel

::

    cp -rf Srivastava None
    cd None
    mv Srivastava.C None.C
    mv Srivastava.H None.H
    
Replace words:

::

    sed -i s/Srivastava/None/g None.H
    sed -i s/Srivastava/None/g None.C

    
Also add to list of diameter models in Make/files    
    
Added this scheme to fvSchemes:

::

    "div(sigmaFDev.particles)"      Gauss limitedLinear 1;   
    

How to scale mesh using OpenFOAM?
---------------------------------

::
    
    transformPoints -scale '(0.001 0.001 0.001)'
    

How to create sets?
-------------------

Define them in the dictionary topoSetDict

::

    actions
    (
        {
            name    extractedFaces;  
            type    faceSet;
            action  new;
            source  boxToFace;
            sourceInfo
            {
                box (-0.00127 0.02855 -0.00127) (0.00127 0.02865 0.00127);   
            }
        }
    
    
        // creation of the faceZone (from the above faceSet)
        {
            name    massFlowSurface;
            type    faceZoneSet;
            action  new;
            source  setAndNormalToFaceZone;
            sourceInfo
            {
                faceSet extractedFaces;
                normal (0    -1    0);
            }
        }
    );

Then run the utility:

::

    topoSet

    
How to view sets in Paraview?
-----------------------------

::

    $ five
    $ paraFoam
    
Turn on Include sets

How to map fields?
------------------

::

    $ cp -rf 14_06_2018_higher_orifice/ 14_06_2018_map_fields
    $ cd 14_06_2018_map_fields/
    $ rm -rf processor*
    $ rm -rf postProcessing
    $ mapFields ../14_06_2018_higher_orifice -parallelSource -sourceTime latestTime
    $ nohup sh run_script_2_continuation.sh > log &
    
    
How to run scripts in series on emps machines?
----------------------------------------------

Create overall run script

::

    # Run script for OpenFOAM 5.0

    cd 026_16mm_55micron_2600kgm3
    nohup sh run_script_2.sh > log &
    echo "Wait until first one is done"
    date
    wait
    date
    echo "First one is done"
    cd ../027_18mm_55micron_2600kgm3
    nohup sh run_script_2.sh > log &
    echo "Second one is done"

Make sure run_script_2.sh is not run in background:

::

    # Run script for OpenFOAM 5.0

    rm -rf log.*

    # Renumber
    renumberMesh -overwrite > log.renumberMesh 2>&1

    # Check the mesh
    checkMesh > log.checkMesh 2>&1

    # Set fields
    rm -rf 0/alpha.particles
    cp 0/alpha.particles.orig 0/alpha.particles
    setFields > log.setFields 2>&1

    # Decompose
    decomposePar > log.decomposePar 2>&1

    # Run
    nohup mpirun -np 6 twoPhaseEulerSrivastavaFoam -parallel > log.twoPhaseEulerSrivastavaFoam 2>&1

Now run the overall script:

::

    nohup sh run.sh > log &



Error in checkMesh due to planes not being aligned
--------------------------------------------------

Number of edges not aligned with or perpendicular to non-empty directions


in cfMesh empty is auto-assigned to all faces - use rename boundary in cfMesh dictionary
Set all tolerances in file>properties to 1e-9 just in case.
Also in pointwise you must project to a database entity, e.g. a cube if the connectors are out of line This with create lines on DB (a domain cannot be created with this), so translate, delete the DB and lines on DB and then translate back again (!)




SnappyHexMesh of Cube
---------------------



Create new solver - twoPhaseEulerSrivastavaSampleFoam
-----------------------------------------------------
 
 
 Copy the code from source

::

    five
    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase
    cp -rf twoPhaseEulerSrivastavaFoam twoPhaseEulerSrivastavaSampleFoam
    cd twoPhaseEulerSrivastavaSampleFoam
    mv twoPhaseEulerSrivastavaFoam.C twoPhaseEulerSrivastavaSampleFoam.C
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g twoPhaseEulerSrivastavaSampleFoam.C
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g Make/files
    ./Allwclean
    ./Allwmake

Added to twoPhaseEulerSampleFoam.C (before int main(): 

::

	#include "cuttingPlane.H"
	#include "sampledCuttingPlane.H"
	#include "sampledPlane.H"
	#include "cellSet.H"

to twoPhaseEulerSampleFoam.C

To /Make/options:

::

    -I$(LIB_SRC)/triSurface/lnInclude\
    -I$(LIB_SRC)/surfMesh/lnInclude \
    -I$(LIB_SRC)/dynamicMesh/lnInclude \


    -ltriSurface\
    -lsurfMesh \
    -ldynamicMesh

Create plane before time loop:

::

    Info<< "\nCreate plane\n" << endl;
    
    point pnt1(-0.00127,0.0216,0);
    point pnt2(0.00127,0.0216,0);
    point pnt3(0.00127,0.0216,0.001);

    plane pl1(pnt1,pnt2,pnt3);
    
    
Create new solver - explicit friction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Solver

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase
    cp -rf twoPhaseEulerSrivastavaFrictionFoam twoPhaseEulerSrivastavaExplicitFrictionFoam
    cd twoPhaseEulerSrivastavaExplicitFrictionFoam
    mv twoPhaseEulerSrivastavaFrictionFoam.C twoPhaseEulerSrivastavaExplicitFrictionFoam.C
    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g twoPhaseEulerSrivastavaExplicitFrictionFoam.C
    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g Make/files
    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g Make/options
    
    
* Interfacial Model               

::

    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g interfacialModels/Make/options
    
    
* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g phaseCompressibleTurbulenceModels/Make/options
           
        
* Phase Model    

::

    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerSrivastavaFrictionFoam/twoPhaseEulerSrivastavaExplicitFrictionFoam/g twoPhaseSystem/Make/options
  

Do the same for twoPhaseEulerSrivastavaFoam
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Do the same for twoPhaseEulerSrivastavaSampleFoam
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Solver

::

    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g Make/options
    
    
* Interfacial Model               

::

    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g interfacialModels/Make/options
    
    
* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g phaseCompressibleTurbulenceModels/Make/options
           
        
* Phase Model    

::

    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerSrivastavaSampleFoam/g twoPhaseSystem/Make/options
    
    

Copy solver from 2.2.2 in order to check gravity term works as it worked for 2.1.0     
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    222
    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase
    cp -rf $FOAM_APP/solvers/multiphase/twoPhaseEulerFoam/ .
    mv twoPhaseEulerFoam twoPhaseEulerGravityVibrationFoam
    cd twoPhaseEulerGravityVibrationFoam
    mv twoPhaseEulerFoam.C twoPhaseEulerGravityVibrationFoam.C
    
* Solver Make files    
   
files   
   
::

    twoPhaseEulerGravityVibrationFoam.C

    EXE = $(FOAM_USER_APPBIN)/twoPhaseEulerGravityVibrationFoam

options    
    
::
    
    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
        -IturbulenceModel \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/kineticTheoryModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude \
        -Iaveraging

    EXE_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -ltwoPhaseEulerGravityVibrationFoamEulerianInterfacialModels \
        -lfiniteVolume \
        -lmeshTools \
        -lincompressibleTransportModels \
        -ltwoPhaseEulerGravityVibrationFoamPhaseModel \
        -ltwoPhaseEulerGravityVibrationFoamKineticTheoryModel
    
    
+ wmake
Making dependency list for source file twoPhaseEulerGravityVibrationFoam.C
SOURCE=twoPhaseEulerGravityVibrationFoam.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/transportModels/incompressible/lnInclude -IturbulenceModel -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/kineticTheoryModels/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -Iaveraging -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/twoPhaseEulerGravityVibrationFoam.o
In file included from twoPhaseEulerGravityVibrationFoam.C:64:0:
/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude/readTimeControls.H: In function ‘int main(int, char**)’:
/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude/readTimeControls.H:38:8: warning: unused variable ‘maxDeltaT’ [-Wunused-variable]
 scalar maxDeltaT =
        ^
g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/transportModels/incompressible/lnInclude -IturbulenceModel -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/kineticTheoryModels/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -Iaveraging -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -Xlinker --add-needed -Xlinker --no-as-needed Make/linux64GccDPOpt/twoPhaseEulerGravityVibrationFoam.o -L/home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib \
     -L/home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib -ltwoPhaseEulerGravityVibrationFoamEulerianInterfacialModels -lfiniteVolume -lmeshTools -lincompressibleTransportModels -ltwoPhaseEulerGravityVibrationFoamPhaseModel -ltwoPhaseEulerGravityVibrationFoamKineticTheoryModel -lOpenFOAM -ldl   -lm -o /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/bin/twoPhaseEulerGravityVibrationFoam
    
* Interfacial Model      

files

::

    dragModels/dragModel/dragModel.C
    dragModels/dragModel/newDragModel.C
    dragModels/Ergun/Ergun.C
    dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C
    dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C
    dragModels/SchillerNaumann/SchillerNaumann.C
    dragModels/Gibilaro/Gibilaro.C
    dragModels/WenYu/WenYu.C
    dragModels/SyamlalOBrien/SyamlalOBrien.C

    LIB = $(FOAM_USER_LIBBIN)/libtwoPhaseEulerGravityVibrationFoamEulerianInterfacialModels

options    
    
::
    
    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -ltwoPhaseEulerGravityVibrationFoamPhaseModel

+ wmake libso interfacialModels
wmakeLnInclude: linking include files to ./lnInclude
Making dependency list for source file dragModels/dragModel/dragModel.C
Making dependency list for source file dragModels/Ergun/Ergun.C
Making dependency list for source file dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C
Making dependency list for source file dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C
Making dependency list for source file dragModels/dragModel/newDragModel.C
Making dependency list for source file dragModels/Gibilaro/Gibilaro.C
Making dependency list for source file dragModels/SchillerNaumann/SchillerNaumann.C
Making dependency list for source file dragModels/WenYu/WenYu.C
Making dependency list for source file dragModels/SyamlalOBrien/SyamlalOBrien.C
SOURCE=dragModels/dragModel/dragModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/dragModel.o
SOURCE=dragModels/dragModel/newDragModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/newDragModel.o
SOURCE=dragModels/Ergun/Ergun.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/Ergun.o
SOURCE=dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/GidaspowErgunWenYu.o
SOURCE=dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/GidaspowSchillerNaumann.o
SOURCE=dragModels/SchillerNaumann/SchillerNaumann.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/SchillerNaumann.o
SOURCE=dragModels/Gibilaro/Gibilaro.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/Gibilaro.o
SOURCE=dragModels/WenYu/WenYu.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/WenYu.o
SOURCE=dragModels/SyamlalOBrien/SyamlalOBrien.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/SyamlalOBrien.o
'/home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib/libtwoPhaseEulerGravityVibrationFoamEulerianInterfacialModels.so' is up to date.
       
        
* Kinetic Theory Model  

files

::

    kineticTheoryModel/kineticTheoryModel.C

    viscosityModel/viscosityModel/viscosityModel.C
    viscosityModel/viscosityModel/newViscosityModel.C
    viscosityModel/Gidaspow/GidaspowViscosity.C
    viscosityModel/Syamlal/SyamlalViscosity.C
    viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C
    viscosityModel/none/noneViscosity.C

    conductivityModel/conductivityModel/conductivityModel.C
    conductivityModel/conductivityModel/newConductivityModel.C
    conductivityModel/Gidaspow/GidaspowConductivity.C
    conductivityModel/Syamlal/SyamlalConductivity.C
    conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C

    radialModel/radialModel/radialModel.C
    radialModel/radialModel/newRadialModel.C
    radialModel/CarnahanStarling/CarnahanStarlingRadial.C
    radialModel/LunSavage/LunSavageRadial.C
    radialModel/SinclairJackson/SinclairJacksonRadial.C

    granularPressureModel/granularPressureModel/granularPressureModel.C
    granularPressureModel/granularPressureModel/newGranularPressureModel.C
    granularPressureModel/Lun/LunPressure.C
    granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C

    frictionalStressModel/frictionalStressModel/frictionalStressModel.C
    frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C
    frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C
    frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C

    LIB = $(FOAM_USER_LIBBIN)/libtwoPhaseEulerGravityVibrationFoamKineticTheoryModel

options    

::

    EXE_INC = \
        -I$(LIB_SRC)/foam/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude

        
result

+ wmake libso kineticTheoryModels
wmakeLnInclude: linking include files to ./lnInclude
Making dependency list for source file kineticTheoryModel/kineticTheoryModel.C
Making dependency list for source file viscosityModel/viscosityModel/viscosityModel.C
Making dependency list for source file viscosityModel/viscosityModel/newViscosityModel.C
Making dependency list for source file viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C
Making dependency list for source file viscosityModel/none/noneViscosity.C
Making dependency list for source file viscosityModel/Syamlal/SyamlalViscosity.C
Making dependency list for source file viscosityModel/Gidaspow/GidaspowViscosity.C
Making dependency list for source file conductivityModel/conductivityModel/conductivityModel.C
Making dependency list for source file conductivityModel/Gidaspow/GidaspowConductivity.C
Making dependency list for source file conductivityModel/conductivityModel/newConductivityModel.C
Making dependency list for source file conductivityModel/Syamlal/SyamlalConductivity.C
Making dependency list for source file conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C
Making dependency list for source file radialModel/radialModel/radialModel.C
Making dependency list for source file radialModel/radialModel/newRadialModel.C
Making dependency list for source file radialModel/CarnahanStarling/CarnahanStarlingRadial.C
Making dependency list for source file radialModel/LunSavage/LunSavageRadial.C
Making dependency list for source file radialModel/SinclairJackson/SinclairJacksonRadial.C
Making dependency list for source file granularPressureModel/granularPressureModel/granularPressureModel.C
Making dependency list for source file granularPressureModel/granularPressureModel/newGranularPressureModel.C
Making dependency list for source file granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C
Making dependency list for source file granularPressureModel/Lun/LunPressure.C
Making dependency list for source file frictionalStressModel/frictionalStressModel/frictionalStressModel.C
Making dependency list for source file frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C
Making dependency list for source file frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C
Making dependency list for source file frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C
SOURCE=kineticTheoryModel/kineticTheoryModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/kineticTheoryModel.o
SOURCE=viscosityModel/viscosityModel/viscosityModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/viscosityModel.o
SOURCE=viscosityModel/viscosityModel/newViscosityModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/newViscosityModel.o
SOURCE=viscosityModel/Gidaspow/GidaspowViscosity.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/GidaspowViscosity.o
SOURCE=viscosityModel/Syamlal/SyamlalViscosity.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/SyamlalViscosity.o
SOURCE=viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/HrenyaSinclairViscosity.o
SOURCE=viscosityModel/none/noneViscosity.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/noneViscosity.o
SOURCE=conductivityModel/conductivityModel/conductivityModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/conductivityModel.o
SOURCE=conductivityModel/conductivityModel/newConductivityModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/newConductivityModel.o
SOURCE=conductivityModel/Gidaspow/GidaspowConductivity.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/GidaspowConductivity.o
SOURCE=conductivityModel/Syamlal/SyamlalConductivity.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/SyamlalConductivity.o
SOURCE=conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/HrenyaSinclairConductivity.o
SOURCE=radialModel/radialModel/radialModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/radialModel.o
SOURCE=radialModel/CarnahanStarling/CarnahanStarlingRadial.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/CarnahanStarlingRadial.o
SOURCE=radialModel/radialModel/newRadialModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/newRadialModel.o
SOURCE=radialModel/LunSavage/LunSavageRadial.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/LunSavageRadial.o
SOURCE=radialModel/SinclairJackson/SinclairJacksonRadial.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/SinclairJacksonRadial.o
SOURCE=granularPressureModel/granularPressureModel/granularPressureModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/granularPressureModel.o
SOURCE=granularPressureModel/granularPressureModel/newGranularPressureModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/newGranularPressureModel.o
SOURCE=granularPressureModel/Lun/LunPressure.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/LunPressure.o
SOURCE=granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/SyamlalRogersOBrienPressure.o
SOURCE=frictionalStressModel/frictionalStressModel/frictionalStressModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/frictionalStressModel.o
SOURCE=frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/newFrictionalStressModel.o
SOURCE=frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/JohnsonJacksonFrictionalStress.o
SOURCE=frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/foam/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/phaseModel/lnInclude -I/home/apr207/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/interfacialModels/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/SchaefferFrictionalStress.o
'/home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib/libtwoPhaseEulerGravityVibrationFoamKineticTheoryModel.so' is up to date.
        
        
* Phase Model 

files

::

    phaseModel/phaseModel.C

    LIB = $(FOAM_USER_LIBBIN)/libtwoPhaseEulerGravityVibrationFoamPhaseModel
    
options    
    
::
    
    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude

    LIB_LIBS = \
        -lincompressibleTransportModels    
    

result    
    
+ wmake libso phaseModel
wmakeLnInclude: linking include files to ./lnInclude
Making dependency list for source file phaseModel/phaseModel.C
SOURCE=phaseModel/phaseModel.C ;  g++ -m64 -Dlinux64 -DWM_DP -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor -O3  -DNoRepository -ftemplate-depth-100 -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/finiteVolume/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/transportModels/incompressible/lnInclude -IlnInclude -I. -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude -I/home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OSspecific/POSIX/lnInclude   -fPIC -c $SOURCE -o Make/linux64GccDPOpt/phaseModel.o
'/home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib/libtwoPhaseEulerGravityVibrationFoamPhaseModel.so' is up to date. 
 
 
Add Gidaspow radial model     
~~~~~~~~~~~~~~~~~~~~~~~~~

::

    cd ~/OpenFOAM/apr207-2.2.2/applications/solvers/multiphase/twoPhaseEulerGravityVibrationFoam/kineticTheoryModels/radialModel
    cp -rf CarnahanStarling/ Gidaspow
    cd Gidaspow
    mv CarnahanStarlingRadial.C GidaspowRadial.C
    mv CarnahanStarlingRadial.H GidaspowRadial.H
    sed -i s/CarnahanStarling/Gidaspow/g GidaspowRadial.C
    sed -i s/CarnahanStarling/Gidaspow/g GidaspowRadial.H
 
 
Edit Gidaspow.C

::

    Foam::tmp<Foam::volScalarField>
    Foam::kineticTheoryModels::radialModels::Gidaspow::g0
    (
        const volScalarField& alpha,
        const dimensionedScalar& alphaMax
    ) const
    {

        return 0.6/(1.0 - pow(alpha/alphaMax, 1.0/3.0));
    }


    Foam::tmp<Foam::volScalarField>
    Foam::kineticTheoryModels::radialModels::Gidaspow::g0prime
    (
        const volScalarField& alpha,
        const dimensionedScalar& alphaMax
    ) const
    {
        return
            (-1.0/5.0)*pow(alpha/alphaMax, -2.0/3.0)
        /(alphaMax*sqr(1.0 - pow(alpha/alphaMax, 1.0/3.0)));
    }
 
* Kinetic Theory Model  

files

::

    kineticTheoryModel/kineticTheoryModel.C

    ...
    radialModel/Gidaspow/GidaspowRadial.C
    ...
 

Effect of mesh size and alignment - hallflow case 
-------------------------------------------------


::

    cp -rf 01_2D_hallflowmeter 10_2D_hallflowmeter
    cd 10_2D_hallflowmeter
    remove_files
    cd constant/polyMesh (remove all files)
    copy all files from 2D_hallflowmeter_unstructured_2 to constant/polyMesh
    transformPoints -scale '(0.001 0.001 0.001)'
    topoSet


    
    
Version 5.0 - new solver - twoPhaseEulerGravityVibrationFoam
------------------------------------------------------------


* Copy folders

::

    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase
    cp -rf twoPhaseEulerSrivastavaFoam twoPhaseEulerGravityVibrationFoam
    cd twoPhaseEulerGravityVibrationFoam
    mv twoPhaseEulerSrivastavaFoam.C twoPhaseEulerGravityVibrationFoam.C
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g twoPhaseEulerGravityVibrationFoam.C


* Solver

::

    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g Make/options
    
    
* Interfacial Model               

::

    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g interfacialModels/Make/options
    
    
* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g phaseCompressibleTurbulenceModels/Make/options
           
        
* Phase Model    

::

    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerSrivastavaFoam/twoPhaseEulerGravityVibrationFoam/g twoPhaseSystem/Make/options

    
OpenFOAM 5.0 - Powder model - axi-symmetric
-------------------------------------------

::

    cp -rf 01_2D_hallflowmeter 12_testing_axisymmetric_grid
    cd 12_testing_axisymmetric_grid
    remove_files
    rm -rf constant/polyMesh/*
    cp -rf /home/apr207/Pointwise_User/2D_hallflowmeter_unstructured_axisymmetric/polyMesh/* ./constant/polyMesh/

    transformPoints -scale '(0.001 0.001 0.001)'

::

    cp -rf 01_2D_hallflowmeter 12_testing_axisymmetric_grid_2
    cd 12_testing_axisymmetric_grid_2
    remove_files
    rm -rf constant/polyMesh/*
    cp -rf /home/apr207/Pointwise_User/2D_hallflowmeter_unstructured_axisymmetric_1_degree/polyMesh/* ./constant/polyMesh/

    transformPoints -scale '(0.001 0.001 0.001)'
    rm -rf 0/*
    cp -rf ../12_testing_axisymmetric_grid/0/* ./0/
    rm -rf system/setFields
    cp -rf ../12_testing_axisymmetric_grid/system/setFieldsDict ./system/
    
    
Version 5.0 - new solver - twoPhaseEulerJohnsonJacksonSchaefferFoam
-------------------------------------------------------------------


* Copy folders

::

    five
    cd $WM_PROJECT_USER_DIR/applications/solvers/multiphase
    cp -rf $FOAM_APP/solvers/multiphase/twoPhaseEulerFoam/ 
    mv twoPhaseEulerFoam twoPhaseEulerJohnsonJacksonSchaefferFoam
    cd twoPhaseEulerJohnsonJacksonSchaefferFoam
    mv twoPhaseEulerFoam.C twoPhaseEulerJohnsonJacksonSchaefferFoam.C
    sed -i s/twoPhaseEulerFoam/twoPhaseEulerJohnsonJacksonSchaefferFoam/g twoPhaseEulerJohnsonJacksonSchaefferFoam.C


* Solver

Make/files

::

    twoPhaseEulerJohnsonJacksonSchaefferFoam.C

    EXE = $(FOAM_USER_APPBIN)/twoPhaseEulerJohnsonJacksonSchaefferFoam
    
Make/options    
      
::

    EXE_INC = \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerJohnsonJacksonSchaefferFoam/phaseCompressibleTurbulenceModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerJohnsonJacksonSchaefferFoam/interfacialModels/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseSystem/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/sampling/lnInclude

    EXE_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie \
        -lturbulenceModels \
        -lcompressibleTurbulenceModels \
        -ltwoPhaseEulerJohnsonJacksonSchaefferFoamPhaseCompressibleTurbulenceModels \
        -lincompressibleTransportModels \
        -ltwoPhaseEulerJohnsonJacksonSchaefferFoamCompressibleTwoPhaseSystem \
        -ltwoPhaseEulerJohnsonJacksonSchaefferFoamCompressibleEulerianInterfacialModels \
        -lfiniteVolume \
        -lfvOptions \
        -lmeshTools \
        -lsampling
     
    
* Interfacial Model               

files

::

    LIB = $(FOAM_USER_LIBBIN)/libtwoPhaseEulerJohnsonJacksonSchaefferFoamCompressibleEulerianInterfacialModels

options    
    
::
    
    EXE_INC = \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/transportModel \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseSystem/lnInclude

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -ltwoPhaseEulerJohnsonJacksonSchaefferFoamCompressibleTwoPhaseSystem \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie

    
    
* Kinetic Theory Model    
     
files     
     
::

    LIB = $(FOAM_USER_LIBBIN)/libtwoPhaseEulerJohnsonJacksonSchaefferFoamPhaseCompressibleFoamTurbulenceModels
    
options

::
    
    EXE_INC = \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerJohnsonJacksonSchaeffer/twoPhaseSystem/lnInclude \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerJohnsonJacksonSchaeffer/interfacialModels/lnInclude\
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/transportModel \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude

    LIB_LIBS = \
        -L$(FOAM_USER_LIBBIN) \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie \
        -lturbulenceModels \
        -lcompressibleTurbulenceModels \
        -lincompressibleTransportModels \
        -ltwoPhaseEulerJohnsonJacksonSchaefferFoamCompressibleTwoPhaseSystem \
        -ltwoPhaseEulerJohnsonJacksonSchaefferFoamCompressibleEulerianInterfacialModels \
        -lfiniteVolume \
        -lfvOptions \
        -lmeshTools
           
        
* Phase Model    

files

::

     LIB = $(FOAM_USER_LIBBIN)/libtwoPhaseEulerJohnsonJacksonSchaefferFoamCompressibleTwoPhaseSystem

options     

::

    EXE_INC = \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseSystem \
        -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerJohnsonJacksonSchaefferFoam/interfacialModels/lnInclude \
        -I$(LIB_SRC)/transportModels/compressible/lnInclude \
        -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
        -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
        -I$(LIB_SRC)/TurbulenceModels/phaseCompressible/lnInclude \
        -I$(LIB_SRC)/finiteVolume/lnInclude \
        -I$(LIB_SRC)/meshTools/lnInclude \
        -I$(LIB_SRC)/sampling/lnInclude

    LIB_LIBS = \
        -lincompressibleTransportModels \
        -lcompressibleTransportModels \
        -lfluidThermophysicalModels \
        -lspecie  
    

Change solver name to twoPhaseEulerDenseFoam
--------------------------------------------

Files

::

    mv twoPhaseEulerJohnsonJacksonSchaefferFoam twoPhaseEulerDenseFoam
    cd twoPhaseEulerDenseFoam
    mv twoPhaseEulerJohnsonJacksonSchaefferFoam.C twoPhaseEulerDenseFoam.C
    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerDenseFoam/g twoPhaseEulerDenseFoam.C


* Solver

::

    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerDenseFoam/g Make/files
    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerDenseFoam/g Make/options
    
    
* Interfacial Model               

::

    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerDenseFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerDenseFoam/g interfacialModels/Make/options
    
    
* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerDenseFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerDenseFoam/g phaseCompressibleTurbulenceModels/Make/options
           
        
* Phase Model    

::

    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerGravityVibrationFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerJohnsonJacksonSchaefferFoam/twoPhaseEulerGravityVibrationFoam/g twoPhaseSystem/Make/options

    
    
    sed -i s/twoPhaseEulerGravityVibrationFoam/twoPhaseEulerDenseFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerGravityVibrationFoam/twoPhaseEulerDenseFoam/g twoPhaseSystem/Make/options


Modify phaseCompressibleTurbulenceModels library to compare Lun and Jenkins-Savage    
    

Change outlet from wall to outlet    
   
::

    36_twoPhaseEulerFoam_srivastava_BCS_JJS_outlet_modified

    transformPoints -scale '(0.001 0.001 0.001)'
    
    topoSet
    
    
Add: 

* JohnsonJacksonSchaefferSrivastava

::

    cp -rf JohnsonJacksonSchaeffer JohnsonJacksonSchaefferSrivastava
    cd JohnsonJacksonSchaefferSrivastava
    mv JohnsonJacksonSchaefferFrictionalStress.C JohnsonJacksonSchaefferSrivastavaFrictionalStress.C
    mv JohnsonJacksonSchaefferFrictionalStress.H JohnsonJacksonSchaefferSrivastavaFrictionalStress.H
    
    sed -i s/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferSrivastava/g JohnsonJacksonSchaefferSrivastavaFrictionalStress.C
    sed -i s/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferSrivastava/g JohnsonJacksonSchaefferSrivastavaFrictionalStress.H


Add selector for wall frictional viscosity in order to avoid oscillations

Export the frictional pressure


    
Case for Morimoto surface
-------------------------

110_denseFoam_amp_0_25mm_len_2_5mm

::

    transformPoints -scale '(0.001 0.001 0.001)'
    


Add new boundary condition JohnsonJacksonParticleThetaVibration
---------------------------------------------------------------

::

    cp -rf JohnsonJacksonParticleTheta JohnsonJacksonParticleThetaVibration
    cd JohnsonJacksonParticleThetaVibration
    mv JohnsonJacksonParticleThetaFvPatchScalarField.C JohnsonJacksonParticleThetaVibrationFvPatchScalarField.C
    mv JohnsonJacksonParticleThetaFvPatchScalarField.H JohnsonJacksonParticleThetaVibrationFvPatchScalarField.H
    
    sed -i s/JohnsonJacksonParticleTheta/JohnsonJacksonParticleThetaVibration/g JohnsonJacksonParticleThetaVibrationFvPatchScalarField.C
    sed -i s/JohnsonJacksonParticleTheta/JohnsonJacksonParticleThetaVibration/g JohnsonJacksonParticleThetaVibrationFvPatchScalarField.H

    
::

    // ADDED: lookup the amplitude
    dimensionedScalar amplitude
    (
        "amplitude",
        dimless,
        db()
       .lookupObject<IOdictionary>
        (
            IOobject::groupName("turbulenceProperties", phased.name())
        )
       .subDict("RAS")
       .subDict("kineticTheoryCoeffs")
       .lookup("amplitude")
    );
    
    
    // ADDED: lookup the frequency
    dimensionedScalar frequency
    (
        "frequency",
        dimless,
        db()
       .lookupObject<IOdictionary>
        (
            IOobject::groupName("turbulenceProperties", phased.name())
        )
       .subDict("RAS")
       .subDict("kineticTheoryCoeffs")
       .lookup("frequency")
    ); 
    


New solver twoPhaseEulerCohesionFoam
------------------------------------

Files

::

    cp -rf twoPhaseEulerDenseFoam twoPhaseEulerCohesionFoam
    cd twoPhaseEulerCohesionFoam
    mv twoPhaseEulerDenseFoam.C twoPhaseEulerCohesionFoam.C
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g twoPhaseEulerCohesionFoam.C


* Solver

::

    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g Make/files
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g Make/options
    
    
* Interfacial Model               

::

    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g interfacialModels/Make/options
    
    
* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g phaseCompressibleTurbulenceModels/Make/options
           
        
* Phase Model    

::

    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerCohesionFoam/g twoPhaseSystem/Make/options

    
Cohesion model
--------------


::

    mkdir cohesionModel
    cp -rf frictionalStressModel/frictionalStressModel cohesionModel
    cp -rf frictionalStressModel/JohnsonJackson cohesionModel
    mv cohesionModel/frictionalStressModel cohesionModel/cohesionModel
    mv cohesionModel/JohnsonJackson cohesionModel/Makkawi

    
::

    cd cohesionModel
    mv frictionalStressModel.C cohesionModel.C  
    mv frictionalStressModel.H cohesionModel.H
    mv newFrictionalStressModel.C newCohesionModel.C
    
::

    sed -i s/frictionalStressModel/cohesionModel/g cohesionModel.C
    sed -i s/frictionalStressModel/cohesionModel/g cohesionModel.H
    sed -i s/frictionalStressModel/cohesionModel/g newCohesionModel.C
    
::

    cd Makkawi
    mv JohnsonJacksonFrictionalStress.C MakkawiCohesion.C  
    mv JohnsonJacksonFrictionalStress.H MakkawiCohesion.H   
    
    
::

    sed -i s/JohnsonJacksonFrictionalStress/MakkawiCohesion/g MakkawiCohesion.C
    sed -i s/JohnsonJacksonFrictionalStress/MakkawiCohesion/g MakkawiCohesion.H
    
    sed -i s/JohnsonJackson/Makkawi/g MakkawiCohesion.C
    sed -i s/JohnsonJackson/Makkawi/g MakkawiCohesion.H
    
    sed -i s/frictionalStressModel/cohesionModel/g MakkawiCohesion.C
    sed -i s/frictionalStressModel/cohesionModel/g MakkawiCohesion.H
    
    

Make/files
    
::

    kineticTheoryModels/cohesionModel/cohesionModel/cohesionModel.C
    kineticTheoryModels/cohesionModel/cohesionModel/newCohesionModel.C
    kineticTheoryModels/cohesionModel/Makkawi/MakkawiCohesion.C


At this point 

::

    ./Allwclean
    ./Allwmake




None cohesion model:


::

    sed -i s/MakkawiCohesion/NoneCohesion/g NoneCohesion.C
    sed -i s/MakkawiCohesion/NoneCohesion/g NoneCohesion.H
    
    sed -i s/Makkawi/None/g NoneCohesion.C
    sed -i s/Makkawi/None/g NoneCohesion.H


    kineticTheoryModels/cohesionModel/None/NoneCohesion.C
    
    
Gidaspow granularPressure model:


::

    sed -i s/LunPressure/GidaspowPressure/g GidaspowPressure.C
    sed -i s/LunPressure/GidaspowPressure/g GidaspowPressure.H
    
    sed -i s/Lun/Gidaspow/g GidaspowPressure.C
    sed -i s/Lun/Gidaspow/g GidaspowPressure.H


    kineticTheoryModels/granularPressureModel/Gidaspow/GidaspowPressure.C

MakkawiSrivastava


::

    sed -i s/MakkawiCohesion/MakkawiSrivastavaCohesion/g MakkawiSrivastavaCohesion.C
    sed -i s/MakkawiCohesion/MakkawiSrivastavaCohesion/g MakkawiSrivastavaCohesion.H
    
    sed -i s/Makkawi/MakkawiSrivastava/g MakkawiSrivastavaCohesion.C
    sed -i s/Makkawi/MakkawiSrivastava/g MakkawiSrivastavaCohesion.H


    kineticTheoryModels/granularPressureModel/Gidaspow/GidaspowPressure.C
    
    
New boundary condition
----------------------


::

    cp -rf JohnsonJacksonParticleTheta JohnsonJacksonParticleCohesion
    mv JohnsonJacksonParticleThetaFvPatchScalarField.C JohnsonJacksonParticleCohesionFvPatchScalarField.C
    mv JohnsonJacksonParticleThetaFvPatchScalarField.H JohnsonJacksonParticleCohesionFvPatchScalarField.H

    sed -i s/JohnsonJacksonParticleTheta/JohnsonJacksonParticleCohesion/g JohnsonJacksonParticleCohesionFvPatchScalarField.C
    sed -i s/JohnsonJacksonParticleTheta/JohnsonJacksonParticleCohesion/g JohnsonJacksonParticleCohesionFvPatchScalarField.H


New solver for sinusoidal gravity oscillation
---------------------------------------------

Files

::

    cp -rf twoPhaseEulerDenseFoam twoPhaseEulerGravityDenseFoam
    cd twoPhaseEulerGravityDenseFoam
    mv twoPhaseEulerDenseFoam.C twoPhaseEulerGravityDenseFoam.C
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g twoPhaseEulerGravityDenseFoam.C

* Solver

::

    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g Make/files
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g Make/options

* Interfacial Model               

::

    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g interfacialModels/Make/options

* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g phaseCompressibleTurbulenceModels/Make/options

* Phase Model    

::

    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerDenseFoam/twoPhaseEulerGravityDenseFoam/g twoPhaseSystem/Make/options





New boundary condition - based on partial slip - change to mixed
----------------------------------------------------------------


::

    cp -rf JohnsonJacksonParticleSlip JohnsonJacksonParticleCohesionMixed
    mv JohnsonJacksonParticleSlipFvPatchVectorField.C JohnsonJacksonParticleCohesionMixedFvPatchVectorField.C
    mv JohnsonJacksonParticleSlipFvPatchVectorField.H JohnsonJacksonParticleCohesionMixedFvPatchVectorField.H

    sed -i s/JohnsonJacksonParticleSlip/JohnsonJacksonParticleCohesionMixed/g JohnsonJacksonParticleCohesionMixedFvPatchVectorField.C
    sed -i s/JohnsonJacksonParticleSlip/JohnsonJacksonParticleCohesionMixed/g JohnsonJacksonParticleCohesionMixedFvPatchVectorField.H

    sed -i s/partialSlipFvPatchVectorField/mixedFvPatchVectorField/g JohnsonJacksonParticleCohesionMixedFvPatchVectorField.C
    sed -i s/partialSlipFvPatchVectorField/mixedFvPatchVectorField/g JohnsonJacksonParticleCohesionMixedFvPatchVectorField.H
    


New solver for centrifugal acceleration
---------------------------------------

Files

::

    cp -rf twoPhaseEulerGravityDenseFoam twoPhaseEulerCentrifugalDenseFoam
    cd twoPhaseEulerCentrifugalDenseFoam
    mv twoPhaseEulerGravityDenseFoam.C twoPhaseEulerCentrifugalDenseFoam.C
    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g twoPhaseEulerCentrifugalDenseFoam.C

* Solver

::

    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g Make/files
    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g Make/options

* Interfacial Model               

::

    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g interfacialModels/Make/options

* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g phaseCompressibleTurbulenceModels/Make/options

* Phase Model    

::

    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerGravityDenseFoam/twoPhaseEulerCentrifugalDenseFoam/g twoPhaseSystem/Make/options




New solver for two-axis rotation
--------------------------------


Files

::

    cp -rf twoPhaseEulerCentrifugalDenseFoam twoPhaseEulerTwoAxisDenseFoam
    cd twoPhaseEulerTwoAxisDenseFoam
    mv twoPhaseEulerCentrifugalDenseFoam.C twoPhaseEulerTwoAxisDenseFoam.C
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g twoPhaseEulerTwoAxisDenseFoam.C

* Solver

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g Make/options

* Interfacial Model               

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g interfacialModels/Make/options

* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g phaseCompressibleTurbulenceModels/Make/options

* Phase Model    

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerTwoAxisDenseFoam/g twoPhaseSystem/Make/options




New solver for one-axis rotation
--------------------------------


Files

::

    cp -rf twoPhaseEulerCentrifugalDenseFoam twoPhaseEulerOneAxisDenseFoam
    cd twoPhaseEulerOneAxisDenseFoam
    mv twoPhaseEulerCentrifugalDenseFoam.C twoPhaseEulerOneAxisDenseFoam.C
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g twoPhaseEulerOneAxisDenseFoam.C

* Solver

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g Make/options

* Interfacial Model               

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g interfacialModels/Make/options

* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g phaseCompressibleTurbulenceModels/Make/options

* Phase Model    

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerOneAxisDenseFoam/g twoPhaseSystem/Make/options




New solver for cohesion and centrifugal force
---------------------------------------------


Files

::

    cp -rf twoPhaseEulerCentrifugalDenseFoam twoPhaseEulerCohesionCentrifugalDenseFoam
    cd twoPhaseEulerCohesionCentrifugalDenseFoam
    mv twoPhaseEulerCentrifugalDenseFoam.C twoPhaseEulerCohesionCentrifugalDenseFoam.C
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g twoPhaseEulerCohesionCentrifugalDenseFoam.C

* Solver

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g Make/options

* Interfacial Model               

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g interfacialModels/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g interfacialModels/Make/options

* Kinetic Theory Model    
     
::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g phaseCompressibleTurbulenceModels/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g phaseCompressibleTurbulenceModels/Make/options

* Phase Model    

::

    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g twoPhaseSystem/Make/files
    sed -i s/twoPhaseEulerCentrifugalDenseFoam/twoPhaseEulerCohesionCentrifugalDenseFoam/g twoPhaseSystem/Make/options








    
    
    
    
    
