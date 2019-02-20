=================
OpenFOAM Language
=================

This is the OpenFOAM language.

.. contents::
   :local:

What are features of C++?
=========================

.. list-table::
   :header-rows: 1
   :widths: 15 60

   * - Feature
     - Meaning
   * - typedefs
     - Alias for a possibly complex type name
   * - function
     - Group of statements that perform a task
   * - pointers
     - Data type that holds addresses to refer to values in memory (e.g. for dynamic memory allocation)    
   * - data structures
     - Data members grouped under one name (e.g. the nodes in a linked list) 
   * - classes
     - Data members and function members grouped under one name
   * - constructor
     - Function member that initialises the instance of it's class
   * - destructor
     - Function member that destroys the instance of it's class    
   * - friends
     - Allows a function or class access to private or protected members of a class
   * - inheritance
     - Allows a class to be created based on another class (so code can be reused)     
   * - virtual member functions
     - Member function that will be redefined in a derived class
   * - abstract class
     - Class that contains at least one virtual function
   * - template
     - Family of classes (class template), functions (function template), variables (variable template) or alias of a family of types (alias template) 
   * - namespace
     - Prevents name conflicts in large projects
     
What is explicit evaluation?
============================

* Evaluate spatial derivatives at the **current** timestep
* Uses **known** values

What is implicit evaluation?
============================
* Evaluate spatial derivatives at **future** timestep
* Uses **unknown** values - generates a matrix equation to be solved          
     
What are the higher level data types?
=====================================

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Type
     - Meaning
   * - ::

           volScalarField

     - scalar, e.g. pressure
   * - ::

           volVectorField

     - vector, e.g. velocity
   * - ::

           volTensorField

     - tensor, e.g. Reynolds Stress    
   * - ::

           surfaceScalarField

     - surface scalar, e.g. flux
   * - ::

           dimensionedScalar

     - constant, e.g. viscosity        
     
     
What are fields?
================

* Arrays of data stored at **cell centres** in the mesh
* Include bouundary information
* Three types 
    - ``volScalarField``
    - ``volVectorField``
    - ``volTensorField``
* Values are stored in named dictionary files in named timestep directories e.g. ``case/0/p`` for pressure
     
What are the three types?
=========================
``<Type>`` refers to:

* ``scalar``
* ``vector``
* ``tensor``
     
What are the five basic classes?
================================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Class
     - Meaning
   * - ::

           fvPatchField
           
       (and derived classes)
           
     - Boundary conditions    
   * - ::

           lduMatrix
           fvMatrix
           
       (and linear solvers)
           
     - Sparse matrices
    
What are the three space-time classes?
======================================

.. list-table::
   :header-rows: 1
   :widths: 10 60 
   
   * - Class
     - Meaning
   * - ::

           polyMesh

     - + Stands for polyhedral mesh
       + Most basic mesh class
       + Contains:
            * ``pointField``
            * ``faceList``
            * ``cellList``
            * ``polyPatchList``
   * - ::

           fvMesh

     - + Extends ``polyMesh`` contains:
            * Cell volumes (``volScalarField``)
            * Cell centres (``volVectorField``)
            * Face area vectors (``surfaceVectorField``)
            * Face centres (``surfaceVectorField``)
            * Face area magnitudes (``surfaceScalarField``)
            * Face motion centres (``surfaceScalarField``)
   * - ::

           Time

     - + Class to control time during OpenFOAM
            * Declared as variable ``runTime``
            * Provides list of saved times ``runTime.times()``
            * Timestep = ``deltaT()``
            * Return current directory name = ``name()``
            * Time increments = ``operator++()`` ``operator+=(scalar)``
            * Write objects to disk = ``write()``
            * Start time, end time = ``startTime()``, ``endTime()``

      
What are the three field algebra classes?
=========================================

.. list-table::
   :header-rows: 1
   :widths: 10 60 
   
   * - Class
     - Meaning
   * - ::

           Field<Type>

     - + Array template class, e.g. ``Field<vector> = vectorField``
       + Renamed using ``typedef`` as:
            * ``scalarField``
            * ``vectorField``
            * ``tensorField``
            * ``symmTensorField``
            * ``tensorThirdField``
            * ``symmTensorThirdField``
   * - ::

           dimensionedField

     - 
   * - ::

           geometricField<Type>

     - + Combination of:
            * ``Field``
            * ``GeometricBoundaryField`` 
            * ``fvMesh``
            * ``dimensionSet``
       + Defines values at all locations in domain with aliases:
            * ``volField<Type>`` 
            * ``surfaceField<Type>``
            * ``pointField<Type>``
     
What are the two discretisation classes?
========================================

.. list-table::
   :header-rows: 1
   :widths: 10 60 
   
   * - Class
     - Meaning
   * - ::

           fvc

     - + Stands for finite volume calculus 
       + Explicit derivative evaluation
       + Input = **known** ``geometricField<Type>``
       + Output = ``geometricField<Type>`` object       
   * - ::

           fvm

     - + Stands for finite volume method
       + Implicit derivative evaluation
       + Input = **unknown** ``geometricField<Type>``
       + Output = ``fvMatrix<Type>`` object, which can be inverted in the matrix equation ``Mx=y``
       

What is a geometricField<Type>?
===============================

* ``volField<Type>``
* ``surfaceField<Type>``
* ``pointField<Type>``

What is the objectRegistry?
===========================

Object registry of entities (dictionaries, fields) which are to be read in or written out

What is the IOobject?
=====================

Defines I/O attributes of entities managed by the object registry.

How is a dictionary object read?
================================

.. list-table::
   :header-rows: 1
   :widths: 30 60 
   
   * - Code
     - Meaning
   * - ::
   
            Info << "Reading transportProperties" << endl;
     - Send message to screen
   * - ::
            
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ
                    IOobject::NO_WRITE
                )
            );
     - Read in at creation
   * - ::
            
            dimensionedScalar nu
            (
                transportProperties.lookup("nu")
            );
     - Lookup viscosity in dictionary

How is volVectorField read from disk?
=====================================

.. list-table::
   :header-rows: 1
   :widths: 30 60 
   
   * - Code
     - Meaning
   * - ::
   
            volVectorField U
            (
                IOobject
                (
                    "U",
                    Times[i].name(),
                    runTime,
                    IOobject::MUST_READ
                ),
                mesh
            )
     - + volVectorField read in from disk
       + Associated with runTime database
       + Must be read

How is volScalarField constructed?
==================================

.. list-table::
   :header-rows: 1
   :widths: 30 60 
   
   * - Code
     - Meaning
   * - ::
   
            volVectorField magU
            (
                IOobject
                (
                    "magU",
                    Times[i].name(),
                    runTime,
                    IOobject::NO_READ
                    IOobject::AUTO_WRITE
                ),
                ::mag(U)
            );
            magU.write();
     - + construct mag(U) object of type volScalarField called magU
       + write it out
       
What are the read and write options?
====================================

.. list-table::
   :header-rows: 1
   :widths: 30 60 
   
   * - Class
     - Meaning
   * - ::
   
            NO_READ
     - Object created
   * - ::
   
            MUST_READ
            READ_IF_PRESENT
     - Object asked to read
   * - ::
   
            NO_WRITE
     - Object destroyed
   * - ::
   
            AUTO_WRITE
            
     - Object asked to write
     
     
How are objects represented in OpenFOAM?
========================================

How to create an object that writes the magnitude of a velocity vector?

.. list-table::
   :header-rows: 1
   :widths: 10 60 
   
   * - Class
     - Meaning
   * - ::

            #include "fvCFD.H"
            int main(int argc, char argv[])
            {
            # include "addTimeOptions.H"
            # include "setRootCase.H"
            # include "createTime.H"
            instantList Times = runTime.times();
            # include "createMesh.H"

     - + One block called main is needed
       + ``#include`` to store commonly used code
       + ``runTime`` is a variable of the OpenFOAM ``Time`` class - for timestepping through code
   * - ::   
   
            for(label i=0; i<runTime.size(); i++)
            {
                runTime.setTime(Times[i],i);
                Info << "Time: " << runTime.value() << endl
                volVectorField U
                (
                    IOobject
                    (
                        "U",
                        Times[i].name(),
                        runTime,
                        IOobject::MUST READ
                    ),
                    mesh
                );
                
     - + Loop over all possible times
       + Read in a ``volVectorField U``
   * - ::   
   
                volScalarField magU
                (
                    IOobject
                    (
                        "magU",
                        Times[i].name(),
                        runTime,
                        IOobject::NO READ,
                        IOobject::AUTO WRITE
                    ),
                    ::mag(U)
                );
                magU.write();
            } return 0;}
            
     - + Construct a ``volScalarField magU``
       + Write out the velocity









       
How is matrix inversion done in fvm?
====================================

* Each operator in fvm constructs particular entries in known ``M`` and ``y`` as a ``fvMatrix`` object
* ``fvMatrix`` is a template class (actual classes are ``fvScalarMatrix`` etc)
* ``fvMatrix`` handles storage via ``lduMatrix`` class
* ``fvMatrix`` class also handles solution

What are lists?
===============

.. list-table::
   :header-rows: 1
   :widths: 10 60 
   
   * - Class
     - Meaning
   * - ::

           List<Type>

     - + Array template class
       + Allows creation of a list of any object of a class e.g. ``List<vector>``
   * - ::

           PtrList<Type>

     - List of pointers
   * - ::

           SLList<Type>

     - Non-intrusive singly-linked list
       
What are fields?
================ 

.. list-table::
   :header-rows: 1
   :widths: 10 60 
   
   * - Class
     - Meaning       
       
   * - ::

           Field<Type>

     - + Array template class, e.g. ``Field<vector> = vectorField``
       + Renamed as ``scalarField, vectorField, tensorField, symmTensorField,``
                    ``tensorThirdField and symmTensorThirdField``
                    
How is memory accessed?
=======================

* Arrays
* Pointers
* References
                    
                    
How is IO Communication done?
=============================

.. list-table::
   :header-rows: 1
   :widths: 10 60 
   
   * - Code
     - Meaning       
       
   * - ::

           Info << "Time = " << runTime.timeName() << nl << endl;

     - Info object is output to the screen

                    
How are derivatives of fields evaluated?
========================================

* Time derivative
* Divergence (div)
    - Spatial derivative
    - Discretised using the flux at the faces
    - e.g. :math:`\nabla \cdot (\mathbf{u}q)` (the advection term)
* Gradient (grad)
    - Spatial derivative
    - e.g. :math:`\nabla p` in the momentum equation
* Laplacian
    - Spatial derivative
    - Discretised as :math:`\nabla \cdot \mu \nabla q`
        + gradient scheme for :math:`\nabla q`
        + interpolation for :math:`\mu`
        + discretisation for :math:`\nabla \cdot`
    - e.g. :math:`\mu \nabla^2 q` in the momentum equation
    
What are the functions for discretisation?
==========================================    
    
.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Meaning
   * - ::

           fvc::ddt(A)
           fvm::ddt(A)
           

     - + Time derivative
       + :math:`\partial A / \partial t`
       + :math:`A` can be a scalar, vector or tensor
   * - ::

           fvc::ddt(rho,A)
           fvm::ddt(rho,A)

     - + Density weighted time derivative
       + :math:`\partial \rho A / \partial t`
       + :math:`\rho` can be any scalar field 
   * - ::

           fvc::d2dt2(rho,A)
           fvm::d2dt2(rho,A)
           
     - + Second density weighted time derivative
       + :math:`\partial / \partial t ( \rho \partial A / \partial t)`   
   * - ::

           fvc::grad(A)
           fvm::grad(A)
           
     - + Gradient
       + :math:`A` can be a scalar or a vector
       + Result is a ``volVectorField`` (from scalar) or a ``volTensorField`` (from vector)
   * - ::

           fvc::div(A)
           fvm::div(A)

     - + Divergence
       + :math:`A` can be a vector or a tensor
       + Result is a ``volScalarField`` (from vector) or a ``volVectorField`` (from tensor) 
   * - ::

           fvc::laplacian(A)
           fvm::laplacian(A)

     - + Laplacian
       + :math:`\nabla^2 A`
   * - ::

           fvc::laplacian(mu, A)
           fvm::laplacian(mu, A)

     - + Laplacian
       + :math:`\nabla \cdot(\mu \nabla A)`    
   * - ::

           fvc::curl(A)
           fvm::curl(A)

     - + Curl
       + :math:`\nabla \times A`        
   * - ::

           fvm::div(phi,A)

     - + Divergence using flux to evaluate this
       + :math:`A` can be a scalar, vector or a tensor   
   * - ::

           fvm::Sp(rho,A)

     - + Implicit evaulation of source term     
   * - ::

           fvm::SuSp(rho,A)

     - + Implicit or explicit evaulation of source term (depending on sign of ``rho``   
     
How can are equations translated into code?
===========================================

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Equation
     - Code
   * - :math:`{{\partial q} \over {\partial t}} + \nabla \cdot q \mathbf{u} = \mu \nabla^2 q`
     - .. code-block:: c

           fvScalarMatrix transport
           (
                fvm::ddt(q)
                + fvm::div(phi,q)
                - fvm::laplacian(mu,q)
           );
           
           // phi is the flux from the momentum equation
           
   * - :math:`{{\partial T} \over {\partial t}} = \kappa \nabla^2 T`
     - .. code-block:: c

           solve(fvm::ddt(T) == kappa*fvc::laplacian(T))
           
           // T is a volScalarField defined on the mesh
           // A discretised representation of the field variable T
           // solve performs matrix inversion for one step
           
   * - :math:`{{\partial k} \over {\partial t}} + \nabla \cdot k \mathbf{u} - \nabla \cdot [(\nu + \nu_t)\nabla k ]=\nu_t[1/2(\nabla \mathbf{u} + \nabla \mathbf{u}^T ]^2 - \varepsilon/k`
     - .. code-block:: c

           solve(
                fvm::ddt(k)
                + fvm::div(phi,k)
                - fvm::laplacian(nu()+nut,k)
                == nut*magSqr(symm(fvc::grad(U)))
                - fvm::Sp(epsilon/k,k)
           );            



How is the PISO algorithm programmed?
=====================================

The PISO (Pressure Implicit with Splitting of Operators) is an efficient method to solve the Navier-Stokes equations

The algorithm can be summed up as follows:

* Set the boundary conditions.
* Solve the discretized momentum equation to compute an intermediate velocity field.
* Compute the mass fluxes at the cells faces.
* Solve the pressure equation.
* Correct the mass fluxes at the cell faces.
* Correct the velocities on the basis of the new pressure field.
* Update the boundary conditions.
* Repeat from 3 for the prescribed number of times.
* Increase the time step and repeat from 1.

The implementation:

* Define the equation for U 

.. code-block:: c

     fvVectorMatrix UEqn
     (
   	fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nu, U)
     );

* Solve the momentum predictor

.. code-block:: c

    solve (UEqn == -fvc::grad(p));
    
* Calculate the ap coefficient and calculate U

.. code-block:: c

    volScalarField rUA = 1.0/UEqn().A();
    U = rUA*UEqn().H();
    
* Calculate the flux

.. code-block:: c

    phi = (fvc::interpolate(U) & mesh.Sf()) 
       + fvc::ddtPhiCorr(rUA, U, phi);
    adjustPhi(phi, U, p);
    
* Define and solve the pressure equation and repeat for the prescribed number of non-orthogonal corrector steps

.. code-block:: c

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUA, p) == fvc::div(phi)
    );
    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();
    
* Correct the flux

.. code-block:: c

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pEqn.flux();
    }
    

* Calculate continuity errors

.. code-block:: c 

    # include "continuityErrs.H"
    

* Perform the momentum corrector step

.. code-block:: c
 
    U -= rUA*fvc::grad(p);
    U.correctBoundaryConditions();
    
    
The following is from the OpenFOAM UK Users Group:
    
What are header files?
======================

Sections of code in separate files that are widely used - all function prototypes in a header file

.. list-table::
   :header-rows: 1
   :widths: 60 60

   * - Equation
     - Code
   * - ::

           #include "CourantNumber.H"
     - File containing code for calculating Courant number      

What is wmake?
==============

* wmake is a make system – directs the compiler to compile specific files in particular ways.
* Controlled through files in Make:
    + files – specifies which user-written files to compile and what to call the result
    + options – specifies additional header files/libraries to be included in the compilation.
  


    
    
    
    
    
    
    
