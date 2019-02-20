twoPhaseEulerFoamErrors
=======================

.. contents::
   :local:

.. highlight:: c++   
   
Error 1 - alpha not showing at t=0
----------------------------------

If setFields does show alpha values at t=0, it is because you need to set the limits of setFields dictionary to exactly correspond with the limits in the model. If you renumberMesh, you must use -overwrite as well, otherwise it will stick the mesh in a funny place 

Adding cells with center within box (-0.04 0 0) (0.04 1 0.01)


Error 2 - floating point exception
----------------------------------

Granular pressure, pa, has gone to infinity (1.8e6) near the wall.

Alpha has also gone above 0.65 

Need to plot everything in gnuplot: residuals, pa, Theta, nua

::

    Courant Number mean: 0.0693311 max: 3.65054
    Max Ur Courant Number = 62.5038
    Time = 0.041

    DILUPBiCG:  Solving for alpha, Initial residual = 0.0210031, Final residual = 4.16181e-11, No Iterations 9
    Dispersed phase volume fraction = 0.571186  Min(alpha) = -0.019327  Max(alpha) = 0.995556
    DILUPBiCG:  Solving for alpha, Initial residual = 0.0107219, Final residual = 3.03307e-11, No Iterations 10
    Dispersed phase volume fraction = 0.571186  Min(alpha) = -0.27625  Max(alpha) = 1.34132
    kinTheory: max(Theta) = 48.117
    kinTheory: min(nua) = 5.03589e-13, max(nua) = 0.0344828
    kinTheory: min(pa) = -7.09087e-10, max(pa) = 3.32578e+07
    [5] #0  Foam::error::printStack(Foam::Ostream&) at ??:?
    [5] #1  Foam::sigFpe::sigHandler(int) at ??:?
    [5] #2   in "/lib/x86_64-linux-gnu/libc.so.6"
    [5] #3  Foam::GAMGSolver::scalingFactor(Foam::Field<double>&, Foam::Field<double> const&, Foam::Field<double> const&, Foam::Field<double> const&) const at ??:?
    [5] #4  Foam::GAMGSolver::scalingFactor(Foam::Field<double>&, Foam::lduMatrix const&, Foam::Field<double>&, Foam::FieldField<Foam::Field, double> const&, Foam::UPtrList<Foam::lduInterfaceField const> const&, Foam::Field<double> const&, unsigned char) const at ??:?
    [5] #5  Foam::GAMGSolver::Vcycle(Foam::PtrList<Foam::lduMatrix::smoother> const&, Foam::Field<double>&, Foam::Field<double> const&, Foam::Field<double>&, Foam::Field<double>&, Foam::Field<double>&, Foam::PtrList<Foam::Field<double> >&, Foam::PtrList<Foam::Field<double> >&, unsigned char) const at ??:?
    [5] #6  Foam::GAMGSolver::solve(Foam::Field<double>&, Foam::Field<double> const&, unsigned char) const at ??:?
    [5] #7  Foam::fvMatrix<double>::solve(Foam::dictionary const&) at ??:?
    [5] #8  
    [5]  at ??:?
    [5] #9  __libc_start_main in "/lib/x86_64-linux-gnu/libc.so.6"
    [5] #10  
    [5]  at ??:?
    [emps-andrew:02679] *** Process received signal ***
    [emps-andrew:02679] Signal: Floating point exception (8)
    [emps-andrew:02679] Signal code:  (-6)
    [emps-andrew:02679] Failing at address: 0x3e800000a77
    [emps-andrew:02679] [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x354b0)[0x7fbe424d84b0]
    [emps-andrew:02679] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x38)[0x7fbe424d8428]
    [emps-andrew:02679] [ 2] /lib/x86_64-linux-gnu/libc.so.6(+0x354b0)[0x7fbe424d84b0]
    [emps-andrew:02679] [ 3] /home/apr207/OpenFOAM/OpenFOAM-2.1.0/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZNK4Foam10GAMGSolver13scalingFactorERNS_5FieldIdEERKS2_S5_S5_+0x64)[0x7fbe435c0934]
    [emps-andrew:02679] [ 4] /home/apr207/OpenFOAM/OpenFOAM-2.1.0/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZNK4Foam10GAMGSolver13scalingFactorERNS_5FieldIdEERKNS_9lduMatrixES3_RKNS_10FieldFieldIS1_dEERKNS_8UPtrListIKNS_17lduInterfaceFieldEEERKS2_h+0x9e)[0x7fbe435c0d0e]
    [emps-andrew:02679] [ 5] /home/apr207/OpenFOAM/OpenFOAM-2.1.0/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZNK4Foam10GAMGSolver6VcycleERKNS_7PtrListINS_9lduMatrix8smootherEEERNS_5FieldIdEERKS8_S9_S9_S9_RNS1_IS8_EESD_h+0x1023)[0x7fbe435c3353]
    [emps-andrew:02679] [ 6] /home/apr207/OpenFOAM/OpenFOAM-2.1.0/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZNK4Foam10GAMGSolver5solveERNS_5FieldIdEERKS2_h+0x38c)[0x7fbe435c42bc]
    [emps-andrew:02679] [ 7] /home/apr207/OpenFOAM/OpenFOAM-2.1.0/platforms/linux64GccDPOpt/lib/libfiniteVolume.so(_ZN4Foam8fvMatrixIdE5solveERKNS_10dictionaryE+0x10b)[0x7fbe44e3bcbb]
    [emps-andrew:02679] [ 8] twoPhaseEulerFoam[0x449c6a]
    [emps-andrew:02679] [ 9] /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf0)[0x7fbe424c3830]
    [emps-andrew:02679] [10] twoPhaseEulerFoam[0x450859]
    [emps-andrew:02679] *** End of error message ***
    --------------------------------------------------------------------------
    mpirun noticed that process rank 5 with PID 2679 on node emps-andrew exited on signal 8 (Floating point exception).
    --------------------------------------------------------------------------

Error 3 - polyMesh folder not present in polymesh (Case Conflict)
-----------------------------------------------------------------

processor7/constant/polymesh (Case Conflict)

rename processor7/constant/polyMesh (Case Conflict) to processor7/constant/polyMesh

::

    --> FOAM FATAL ERROR: 
    Cannot find file "points" in directory "polyMesh" in times 0 down to constant

        From function Time::findInstance(const fileName&, const word&, const IOobject::readOption, const word&)
        in file db/Time/findInstance.C at line 203.

Error 4 - cannot find phaseProperties
-------------------------------------

Because using the alias paraview4 invokes openfoam 2.4.0, so to run, you must use 2.2.2 again

::

    [0] --> FOAM FATAL IO ERROR: 
    [0] cannot open file
    [0] 
    [0] file: /home/apr207/OpenFOAM/apr207-2.1.0/run/eulerian_002_bin_discharge/processor0/constant/phaseProperties at line 0.


Error 5 - cannot find file
--------------------------

The sets folder is written when checkMesh employed

::

    [1] --> FOAM FATAL IO ERROR: 
    [1] cannot find file
    [1] 
    [1] file: /home/apr207/OpenFOAM/apr207-2.2.2/run/eulerian_011_analytical_equation/processor1/0/U1 at line 0.

Error 5 - zig-zag results that are bizzare
------------------------------------------

Remove mean files from 0 folder

If there are no mean results - re-run 


Error 6 - no rule
-----------------

::

    ./Allwclean
    ./Allwmake

Make sure that the files are defintely in the correct folder! If not, there will be an error like this:

::

    make: *** No rule to make target 'derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.dep', needed by 'Make/linux64GccDPOpt/dependencies'. Stop.

Error 7 - twoPhaseSystem.H: No such file or directory
-----------------------------------------------------

This error occurs because of renaming of variables and the way they are looked up in the dictinaries between 2.3.x and 2.2.2.

::

    ./Allwclean
    ./Allwmake

::
    
    derivedFvPatchFields/JohnsonJacksonParticleSlip/JohnsonJacksonParticleSlipFvPatchVectorField.C:24:28: fatal error: twoPhaseSystem.H: No such file or directory

Replace header files - twoPhaseSystem.H is used in 2.3.x to allow different names for the two phases, e.g. air and water, looked up in the phaseProperties dictionary.
 
The equivalent in 2.2.2 is phaseModel.H, which looks up everything from transportProperties.

 
.. list-table::
    :header-rows: 1
    :widths: 60 60

    * - Old
      - New
    * - ::

            #include "twoPhaseSystem.H"
            
      - ::

            #include "phaseModel.H"

However, this requires that we change the rest of the code, in order to lookup everything from the new dictionary.           
           
           
.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Old
     - New           
   * - ::
       
            void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::updateCoeffs()
            {
                if (updated())
                {
                    return;
                }

                // lookup the fluid model and the phase
                const twoPhaseSystem& fluid = db().lookupObject<twoPhaseSystem>
                (
                    "phaseProperties"
                );

                const phaseModel& phased
                (
                    fluid.phase1().name() == dimensionedInternalField().group()
                ? fluid.phase1()
                : fluid.phase2()
                );

                // lookup all the fields on this patch
                const fvPatchScalarField& alpha
                (
                    patch().lookupPatchField<volScalarField, scalar>
                    (
                        phased.volScalarField::name()
                    )
                );

                const fvPatchScalarField& gs0
                (
                    patch().lookupPatchField<volScalarField, scalar>
                    (
                        IOobject::groupName("gs0", phased.name())
                    )
                );

                const scalarField nu
                (
                    patch().lookupPatchField<volScalarField, scalar>
                    (
                        IOobject::groupName("nut", phased.name())
                    )
                );

                word ThetaName(IOobject::groupName("Theta", phased.name()));

                const fvPatchScalarField& Theta
                (
                    db().foundObject<volScalarField>(ThetaName)
                ? patch().lookupPatchField<volScalarField, scalar>(ThetaName)
                : alpha
                );

                // lookup the packed volume fraction
                dimensionedScalar alphaMax
                (
                    "alphaMax",
                    dimless,
                    db()
                .lookupObject<IOdictionary>
                    (
                        IOobject::groupName("turbulenceProperties", phased.name())
                    )
                .subDict("RAS")
                .subDict("kineticTheoryCoeffs")
                .lookup("alphaMax")
                );

                // calculate the slip value fraction
                scalarField c
                (
                    constant::mathematical::pi
                *alpha
                *gs0
                *specularityCoefficient_.value()
                *sqrt(3.0*Theta)
                /max(6.0*nu*alphaMax.value(), SMALL)
                );

                this->valueFraction() = c/(c + patch().deltaCoeffs());

                partialSlipFvPatchVectorField::updateCoeffs();
            }
            
     - ::
      
            void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::updateCoeffs()
            {
                if (updated())
                {
                    return;
                }


                // lookup all the fields on this patch
                const fvPatchScalarField& alpha =
                    patch().lookupPatchField<volScalarField, scalar>("alpha1");
                const fvPatchScalarField& nu =
                    patch().lookupPatchField<volScalarField, scalar>("nuEff1");
                const fvPatchScalarField& gs0 =
                    patch().lookupPatchField<volScalarField, scalar>("gs0");
                const fvPatchScalarField& Theta =
                    patch().lookupPatchField<volScalarField, scalar>("Theta");


                const dictionary& kineticTheoryProperties =
                    db()
                .lookupObject<IOdictionary>
                    (
                        "kineticTheoryProperties"
                    );

                const dimensionedScalar alphaMax(kineticTheoryProperties.lookup("alphaMax"));	

                // calculate the slip value fraction
                scalarField c
                (
                    constant::mathematical::pi
                *alpha
                *gs0
                *specularityCoefficient_.value()
                *sqrt(3.0*Theta)
                /max(6.0*nu*alphaMax.value(), SMALL)
                );

                this->valueFraction() = c/(c + patch().deltaCoeffs());

                partialSlipFvPatchVectorField::updateCoeffs();
            }
           
           
Error 8 - can't find patch type
-------------------------------  
   
::

    [2] --> FOAM FATAL IO ERROR:
    [2] Unknown patchField type JohnsonJacksonParticleSlip for patch type wall
  
Add this to controlDict - when you have custom libraries, OpenFOAM needs to be given the name of the new library. However, if the solver has been compiled correctly, this should now be neccessary./

::

   libs ("libkineticTheoryModelSrivastava.so");
   
   
   
Error 9 - duplicate entries
---------------------------
   
::

    Duplicate entry Gidaspow in runtime selection table viscosityModel
    Duplicate entry Syamlal in runtime selection table viscosityModel
    Duplicate entry HrenyaSinclair in runtime selection table viscosityModel
    Duplicate entry none in runtime selection table viscosityModel
    Duplicate entry Gidaspow in runtime selection table conductivityModel
    Duplicate entry CarnahanStarling in runtime selection table radialModel
    Duplicate entry LunSavage in runtime selection table radialModel
    Duplicate entry SinclairJackson in runtime selection table radialModel   
    Duplicate entry Lun in runtime selection table granularPressureModel   
    Duplicate entry SyamlalRogersOBrien in runtime selection table granularPressureModel   
    Duplicate entry JohnsonJackson in runtime selection table frictionalStressModel   
    Duplicate entry Schaeffer in runtime selection table frictionalStressModel   

This is because these two libraries duplicate the entries:    

::

    /home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib/libkineticTheoryModel.so
    /home/apr207/OpenFOAM/apr207-2.2.2/platforms/linux64GccDPOpt/lib/libkineticTheoryModelSrivastava.so   
   
This is because you are running twoPhaseEulerFoam and calling the library twoPhaseEulerSrivastava in controlDict, so these two conflict!
   
   
Error 10 - mpi failure at end
-----------------------------   

This was because you are running twoPhaseEulerFoam and not twoPhaseeulerSrivastava
   
::

    fieldAverage fieldAverage1 output:
        Calculating averages
        Writing average fields

    End

    Finalising parallel run
    [emps-andrew:00860] *** Process received signal ***
    [emps-andrew:00860] Signal: Aborted (6)
    [emps-andrew:00860] Signal code:  (-6)
    [emps-andrew:00860] [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x354b0)[0x7f872ddda4b0]
    [emps-andrew:00860] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x38)[0x7f872ddda428]
    [emps-andrew:00860] [ 2] /lib/x86_64-linux-gnu/libc.so.6(abort+0x16a)[0x7f872dddc02a]
    [emps-andrew:00860] [ 3] /lib/x86_64-linux-gnu/libc.so.6(+0x777ea)[0x7f872de1c7ea]
    [emps-andrew:00860] [ 4] /lib/x86_64-linux-gnu/libc.so.6(+0x80baf)[0x7f872de25baf]
    [emps-andrew:00860] [ 5] /lib/x86_64-linux-gnu/libc.so.6(cfree+0x4c)[0x7f872de2953c]
    [emps-andrew:00860] [ 6] [emps-andrew:00861] *** Process received signal ***
    [emps-andrew:00861] Signal: Aborted (6)
    [emps-andrew:00861] Signal code:  (-6)
    [emps-andrew:00861] [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x354b0)[0x7f181f02b4b0]
    [emps-andrew:00861] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x38)[0x7f181f02b428]
    [emps-andrew:00861] [ 2] /lib/x86_64-linux-gnu/libc.so.6(abort+0x16a)[0x7f181f02d02a]
    /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD1Ev+0xdd)[0x7f872eec55ad]
    [emps-andrew:00860] [ 7] [emps-andrew:00863] *** Process received signal ***
    [emps-andrew:00863] Signal: Aborted (6)
    [emps-andrew:00863] Signal code:  (-6)
    [emps-andrew:00863] [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x354b0)[0x7efc67b424b0]
    [emps-andrew:00863] [ 1] [emps-andrew:00861] [ 3] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x38)[0x7efc67b42428]
    /lib/x86_64-linux-gnu/libc.so.6(+0x777ea)[0x7f181f06d7ea]
    [emps-andrew:00861] [ 4] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam15dictionaryEntryD0Ev+0x2a)[0x7f872eeda88a]
    [emps-andrew:00860] [ 8] [emps-andrew:00863] [ 2] /lib/x86_64-linux-gnu/libc.so.6(+0x80c71)[0x7f181f076c71]
    /lib/x86_64-linux-gnu/libc.so.6(abort+0x16a)[0x7efc67b4402a]
    [emps-andrew:00863] [ 3] [emps-andrew:00861] /lib/x86_64-linux-gnu/libc.so.6(+0x777ea)[0x7efc67b847ea]
    [emps-andrew:00863] [ 4] [ 5] /lib/x86_64-linux-gnu/libc.so.6(+0x80c71)[0x7efc67b8dc71]
    [emps-andrew:00863] [ 5] /lib/x86_64-linux-gnu/libc.so.6(cfree+0x4c)[0x7efc67b9153c]
    /lib/x86_64-linux-gnu/libc.so.6(cfree+0x4c)[0x7f181f07a53c]
    [emps-andrew:00861] [ 6] /lib/x86_64-linux-gnu/libc.so.6(+0x3a025)[0x7f181f030025]
    [emps-andrew:00861] [ 7] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD1Ev+0xdd)[0x7f872eec55ad]
    [emps-andrew:00860] [ 9] /lib/x86_64-linux-gnu/libc.so.6(+0x3a045)[0x7f181f030045]
    [emps-andrew:00861] [ 8] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/openmpi-system/libPstream.so(+0x476a)[0x7f181edef76a]
    [emps-andrew:00861] [ 9] twoPhaseEulerFoam[0x44050a]
    [emps-andrew:00861] [10] [emps-andrew:00863] /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf0)[0x7f181f016830]
    [emps-andrew:00861] [11] twoPhaseEulerFoam[0x441eb9]
    [emps-andrew:00861] *** End of error message ***
    [ 6] /lib/x86_64-linux-gnu/libc.so.6(+0x3a025)[0x7efc67b47025]
    /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD0Ev+0x9)[0x7f872eec5719]
    [emps-andrew:00860] [10] [emps-andrew:00863] [ 7] /lib/x86_64-linux-gnu/libc.so.6(+0x3a045)[0x7efc67b47045]
    [emps-andrew:00863] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam5debug20deleteControlDictPtrD1Ev+0x2b4)[0x7f872ee1b134]
    [emps-andrew:00860] [11] [ 8] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/openmpi-system/libPstream.so(+0x476a)[0x7efc6790676a]
    [emps-andrew:00863] [ 9] /lib/x86_64-linux-gnu/libc.so.6(__cxa_finalize+0x9a)[0x7f872dddf36a]
    [emps-andrew:00860] [12] twoPhaseEulerFoam[0x44050a]
    [emps-andrew:00863] [10] /lib/x86_64-linux-gnu/libc.so.6(__libc_start_main+0xf0)[0x7efc67b2d830]
    [emps-andrew:00863] [11] twoPhaseEulerFoam[0x441eb9]
    [emps-andrew:00863] *** End of error message ***
    /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(+0x1fd623)[0x7f872ee11623]
    [emps-andrew:00860] *** End of error message ***
    [emps-andrew:00864] *** Process received signal ***
    [emps-andrew:00864] Signal: Aborted (6)
    [emps-andrew:00864] Signal code:  (-6)
    [emps-andrew:00864] [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x354b0)[0x7f12485454b0]
    [emps-andrew:00864] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x38)[0x7f1248545428]
    [emps-andrew:00864] [ 2] /lib/x86_64-linux-gnu/libc.so.6(abort+0x16a)[0x7f124854702a]
    [emps-andrew:00864] [ 3] /lib/x86_64-linux-gnu/libc.so.6(+0x777ea)[0x7f12485877ea]
    [emps-andrew:00864] [ 4] /lib/x86_64-linux-gnu/libc.so.6(+0x80c71)[0x7f1248590c71]
    [emps-andrew:00864] [ 5] /lib/x86_64-linux-gnu/libc.so.6(cfree+0x4c)[0x7f124859453c]
    [emps-andrew:00864] [ 6] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD1Ev+0xdd)[0x7f12496305ad]
    [emps-andrew:00864] [ 7] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam15dictionaryEntryD0Ev+0x2a)[0x7f124964588a]
    [emps-andrew:00864] [ 8] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD1Ev+0xdd)[0x7f12496305ad]
    [emps-andrew:00864] [ 9] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD0Ev+0x9)[0x7f1249630719]
    [emps-andrew:00864] [10] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam5debug20deleteControlDictPtrD1Ev+0x2b4)[0x7f1249586134]
    [emps-andrew:00864] [11] /lib/x86_64-linux-gnu/libc.so.6(__cxa_finalize+0x9a)[0x7f124854a36a]
    [emps-andrew:00864] [12] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(+0x1fd623)[0x7f124957c623]
    [emps-andrew:00864] *** End of error message ***
    [emps-andrew:00865] *** Process received signal ***
    [emps-andrew:00865] Signal: Aborted (6)
    [emps-andrew:00865] Signal code:  (-6)
    [emps-andrew:00865] [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x354b0)[0x7fdbfd0364b0]
    [emps-andrew:00865] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x38)[0x7fdbfd036428]
    [emps-andrew:00865] [ 2] /lib/x86_64-linux-gnu/libc.so.6(abort+0x16a)[0x7fdbfd03802a]
    [emps-andrew:00865] [ 3] /lib/x86_64-linux-gnu/libc.so.6(+0x777ea)[0x7fdbfd0787ea]
    [emps-andrew:00865] [ 4] /lib/x86_64-linux-gnu/libc.so.6(+0x7e6ed)[0x7fdbfd07f6ed]
    [emps-andrew:00865] [ 5] /lib/x86_64-linux-gnu/libc.so.6(+0x80678)[0x7fdbfd081678]
    [emps-andrew:00865] [ 6] /lib/x86_64-linux-gnu/libc.so.6(cfree+0x4c)[0x7fdbfd08553c]
    [emps-andrew:00865] [ 7] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD1Ev+0xdd)[0x7fdbfe1215ad]
    [emps-andrew:00865] [ 8] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam15dictionaryEntryD0Ev+0x2a)[0x7fdbfe13688a]
    [emps-andrew:00865] [ 9] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD1Ev+0xdd)[0x7fdbfe1215ad]
    [emps-andrew:00865] [10] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD0Ev+0x9)[0x7fdbfe121719]
    [emps-andrew:00865] [11] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam5debug20deleteControlDictPtrD1Ev+0x2b4)[0x7fdbfe077134]
    [emps-andrew:00865] [12] /lib/x86_64-linux-gnu/libc.so.6(__cxa_finalize+0x9a)[0x7fdbfd03b36a]
    [emps-andrew:00865] [13] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(+0x1fd623)[0x7fdbfe06d623]
    [emps-andrew:00865] *** End of error message ***
    [emps-andrew:00862] *** Process received signal ***
    [emps-andrew:00862] Signal: Aborted (6)
    [emps-andrew:00862] Signal code:  (-6)
    [emps-andrew:00862] [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x354b0)[0x7f55f39fb4b0]
    [emps-andrew:00862] [ 1] /lib/x86_64-linux-gnu/libc.so.6(gsignal+0x38)[0x7f55f39fb428]
    [emps-andrew:00862] [ 2] /lib/x86_64-linux-gnu/libc.so.6(abort+0x16a)[0x7f55f39fd02a]
    [emps-andrew:00862] [ 3] /lib/x86_64-linux-gnu/libc.so.6(+0x777ea)[0x7f55f3a3d7ea]
    [emps-andrew:00862] [ 4] /lib/x86_64-linux-gnu/libc.so.6(+0x80c71)[0x7f55f3a46c71]
    [emps-andrew:00862] [ 5] /lib/x86_64-linux-gnu/libc.so.6(cfree+0x4c)[0x7f55f3a4a53c]
    [emps-andrew:00862] [ 6] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD1Ev+0xdd)[0x7f55f4ae65ad]
    [emps-andrew:00862] [ 7] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam15dictionaryEntryD0Ev+0x2a)[0x7f55f4afb88a]
    [emps-andrew:00862] [ 8] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD1Ev+0xdd)[0x7f55f4ae65ad]
    [emps-andrew:00862] [ 9] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam10dictionaryD0Ev+0x9)[0x7f55f4ae6719]
    [emps-andrew:00862] [10] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(_ZN4Foam5debug20deleteControlDictPtrD1Ev+0x2b4)[0x7f55f4a3c134]
    [emps-andrew:00862] [11] /lib/x86_64-linux-gnu/libc.so.6(__cxa_finalize+0x9a)[0x7f55f3a0036a]
    [emps-andrew:00862] [12] /home/apr207/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/lib/libOpenFOAM.so(+0x1fd623)[0x7f55f4a32623]
    [emps-andrew:00862] *** End of error message ***
    --------------------------------------------------------------------------
    mpirun noticed that process rank 2 with PID 862 on node emps-andrew exited on signal 6 (Aborted).
    --------------------------------------------------------------------------
   
   
Error 11 - member function not implemented
------------------------------------------   

This was because you are running twoPhaseEulerFoam and not twoPhaseEulerSrivastava
   
   
   
 ::

    In file included from /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/token.H:49:0,
                    from /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/UILListIO.C:28,
                    from /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/UILList.C:92,
                    from /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/UILList.H:330,
                    from /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/ILList.H:39,
                    from /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/IDLList.H:35,
                    from /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/entry.H:45,
                    from /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/dictionary.H:53,
                    from lnInclude/frictionalStressModel.H:35,
                    from frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.H:36,
                    from frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.C:26:
    lnInclude/frictionalStressModel.H: In instantiation of ‘static Foam::autoPtr<Foam::frictionalStressModel> Foam::frictionalStressModel::adddictionaryConstructorToTable<frictionalStressModelType>::New(const Foam::dictionary&) [with frictionalStressModelType = Foam::JohnsonJacksonSchaefferFrictionalStress]’:
    lnInclude/frictionalStressModel.H:73:5:   required from ‘Foam::frictionalStressModel::adddictionaryConstructorToTable<frictionalStressModelType>::adddictionaryConstructorToTable(const Foam::word&) [with frictionalStressModelType = Foam::JohnsonJacksonSchaefferFrictionalStress]’
    frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.C:35:5:   required from here
    /home/apr207/OpenFOAM/OpenFOAM-2.2.2/src/OpenFOAM/lnInclude/runTimeSelectionTables.H:76:66: error: invalid new-expression of abstract class type ‘Foam::JohnsonJacksonSchaefferFrictionalStress’
                return autoPtr< baseType >(new baseType##Type parList);           \
                                                                    ^
    lnInclude/frictionalStressModel.H:73:5: note: in expansion of macro ‘declareRunTimeSelectionTable’
        declareRunTimeSelectionTable
     ^

::

    In file included from frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.C:26:0:
    frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.H:47:7: note:   because the following virtual functions are pure within ‘Foam::JohnsonJacksonSchaefferFrictionalStress’:
    class JohnsonJacksonSchaefferFrictionalStress
        ^
        
::

    In file included from frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.H:36:0,
                    from frictionalStressModel/JohnsonJacksonSchaeffer/JohnsonJacksonSchaefferFrictionalStress.C:26:
    lnInclude/frictionalStressModel.H:125:37: note: 	virtual Foam::tmp<Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> > Foam::frictionalStressModel::muf(const volScalarField&, const dimensionedScalar&, const volScalarField&, const volSymmTensorField&, const dimensionedScalar&) const
            virtual tmp<volScalarField> muf
                                        ^

There is nothing unclear about the error message. Your class Foam::JohnsonJacksonSchaefferFrictionalStress has at least one member that is not implemented, which means it is abstract. You cannot instantiate an abstract class.

Where to look: frictionModel.H needs a new defintion for muf (and so do all the other friction models):

::

        virtual tmp<volScalarField> muf
        (
            const volScalarField& alpha1,
            const volScalarField& Theta,
            const dimensionedScalar& alphaMinFriction,
            const dimensionedScalar& alphaMax,
            const volScalarField& pf,
            const volSymmTensorField& D,
            const dimensionedScalar& phi
        ) const = 0;



Error 12 - Not enough arguments
-------------------------------

::

    In member function ‘void Foam::kineticTheoryModel::solve(const volTensorField&)’:
    kineticTheoryModel/kineticTheoryModel.C:369:9: error: no matching function for call to ‘Foam::frictionalStressModel::muf(const volScalarField&, const dimensionedScalar&, Foam::volScalarField&, Foam::volSymmTensorField&, const dimensionedScalar&)’
            )
            ^
    In file included from kineticTheoryModel/kineticTheoryModel.H:44:0,
                    from kineticTheoryModel/kineticTheoryModel.C:26:
    lnInclude/frictionalStressModel.H:125:37: note: candidate: virtual Foam::tmp<Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> > Foam::frictionalStressModel::muf(const volScalarField&, const volScalarField&, const dimensionedScalar&, const dimensionedScalar&, const volScalarField&, const volSymmTensorField&, const dimensionedScalar&) const
            virtual tmp<volScalarField> muf
                                        ^
    lnInclude/frictionalStressModel.H:125:37: note:   candidate expects 7 arguments, 5 provided            
            

Simply add more arguments to match the .H files:            
            
::

    // frictional shear stress, Eq. 3.30, p. 52
    volScalarField muf
    (
        frictionalStressModel_->muf
        (
            alpha1_,
            Theta_,
            alphaMinFriction_,
            alphaMax_,
            pf,
            D,
            phi_
        )
    );  
   
   
Error 13 - OpenFOAM version 5.0 - cannot open case directory 
------------------------------------------------------------

blockMesh has not been run

::
   
    [0] 
    [0] 
    [0] --> FOAM FATAL ERROR: 
    [0] twoPhaseEulerFoam: cannot open case directory "/home/apr207/OpenFOAM/apr207-5.0/run/fluidisedBed/processor0"
    [0] 
    [0] 
    FOAM parallel run exiting
    [0]  
   

Error 14 - crazy zig-zag results
--------------------------------

Parallization probably not working properly - changed run script to use run functions

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

    echo
    echo "Creating files for paraview post-processing"
    echo
    touch foam.foam
   
   
Error 15 - negative initial temperature
---------------------------------------   

::
   
    [0] --> FOAM FATAL ERROR: 
    [0] Negative initial temperature T0: -730.68
    [0]    
   
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

Error 16 - cannot find file
---------------------------

Remove mean files

::

    --> FOAM FATAL IO ERROR:
    [2] cannot find file
    [2]
    [2] file: /home/apr207/OpenFOAM/apr207-2.2.2/run/eulerian_043_walls_johnson_jackson_schaeffer/processor2/0/U1 at line 0.
    [2]
    [2]     From function regIOobject::readStream()
    [2]     in file db/regIOobject/regIOobjectRead.C at line 73.
  

  
Error 17 - Simulation stopped after I logged out
------------------------------------------------
  
::

    nohup sh run_script.sh &

Or to direct it to a log file:

::

    nohup sh run_script.sh > log &
  
  
::

    # Run script for OpenFOAM 5.0

    rm -rf log.*

    # Check the mesh
    checkMesh

    # Set fields
    rm -rf 0/alpha.particles
    cp 0/alpha.particles.orig 0/alpha.particles
    setFields

    # Decompose
    decomposePar

    # Renumber
    renumberMesh -overwrite

    # Run
    nohup nice -19 mpirun -np 6 twoPhaseEulerFoam -parallel > log 2>&1

    # Reconstruct
    reconstructPar

    rm -r processor*

    echo "Creating files for paraview post-processing"
    touch foam.foam
  
  
  
Error 18 - Schaeffer model crashes
----------------------------------

Must change alphaMinFriction so that alphaMax - alphaMinFriction = around 0.02

So if alphaMax = 0.65, then alphaMinFriction = 0.63

::

    ---
    #0 Foam::error::printStack(Foam::Ostream&) at ??:?
    #1 Foam::sigFpe::sigHandler(int) at ??:?
    #2 in "/lib/x86_64-linux-gnu/libc.so.6"
    #3 in "/lib/x86_64-linux-gnu/libm.so.6"
    #4 pow in "/lib/x86_64-linux-gnu/libm.so.6"
    #5 Foam::pow(Foam::Field<double>&, Foam::UList<double> const&, double const&) at ??:?
    #6 Foam::tmp<Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> > Foam::pow<Foam::fvPatchField, Foam::volMesh>(Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> const&, Foam::dimensioned<double> const&) at ??:?
    #7 Foam::dragModels::WenYu::CdRe() const at ??:?
    #8 Foam::dragModels::GidaspowErgunWenYu::CdRe() const at ??:?
    #9 Foam::dragModel::K() const at ??:?
    #10 Foam::BlendedInterfacialModel<Foam::dragModel>::K() const at ??:?
    #11 Foam::twoPhaseSystem::dragCoeff() const at ??:?
    #12 
    at ??:?
    #13 __libc_start_main in "/lib/x86_64-linux-gnu/libc.so.6"
    #14 
    at ??:?
    ---


The Schaeffer model generates an exponentially increasing frictional pressure for particulate phase fractions above the specified alphaMinFriction of 0.5. The RAS setup also includes a particle pressure which allows the particulate fraction to reach 0.62. At this value, the frictional pressure generated by the Schaeffer model is very large, and presumably not physical. You've already found the answer to the instability; increasing alphaMinFriction.

Models work with some coefficients, and are unstable with others. If you have data which suggests the Schaeffer model is wrong, or you can spot something incorrect in the implementation, then please feel free to re-open this report. As that information is currently lacking, I'm marking this as closed.

NB: Re the CFD Online thread, frictionalPressurePrime is the derivative of the frictional pressure with respect to the volume fraction. The Schaeffer model is implemented consistently in this respect.




Error 18 - Simulation stopped after I closed window
---------------------------------------------------

Implies it was running in the foreground: add & to end of nohup to run in background and get rid of nice

::
  
    # Run script for OpenFOAM 5.0

    rm -rf log.*

    # Check the mesh
    checkMesh

    # Set fields
    rm -rf 0/alpha.particles
    cp 0/alpha.particles.orig 0/alpha.particles
    setFields

    # Decompose
    decomposePar

    # Renumber
    renumberMesh -overwrite

    # Run
    nohup mpirun -np 6 twoPhaseEulerFoam -parallel > log 2>&1 &

    # Reconstruct
    #reconstructPar

    #rm -r processor*

    #echo "Creating files for paraview post-processing"
    #touch foam.foam

    
    
Error 19 - cannot find p file
-----------------------------

There is a syntax error in the BCs description, e.g. missing semicolon

::

    walls
    {
        type            fixedValue
        value           uniform (0 0 0);
    }


::

    [0] --> FOAM FATAL ERROR:
    [0] cannot find file "/home/apr207/OpenFOAM/apr207-5.0/run/020_srivastava_syamlal_jjs_friction_inlet_outlet_walls/processor0/0/p"
    [0]
    [0]     From function virtual Foam::autoPtr<Foam::ISstream> Foam::fileOperations::uncollatedFileOperation::readStream(Foam::regIOobject&, con$
    [0]     in file global/fileOperations/uncollatedFileOperation/uncollatedFileOperation.C at line 505.


    
Error 20 - cannot find case directory
-------------------------------------

There was no domain decomposition because of an error in decomposePar
    
::
    
    [0] --> FOAM FATAL ERROR:
    [0] twoPhaseEulerFoam: cannot open case directory "/home/apr207/OpenFOAM/apr207-5.0/run/025_srivastava_inletOutlet_prghPressure_pressureInletOutletVelocity_prghTotalPressure/processor0"

Error 21 - Theta and da not available, look them up
---------------------------------------------------

::

    volScalarField Theta = pf.mesh().lookupObject<volScalarField>("Theta");
    volScalarField da = pf.mesh().lookupObject<volScalarField>("da");    
    
    
Error 22 - invalid types
------------------------

Its because you used this:

::

    D.boundaryField[patchi]

Instead: 

::

    D.boundaryField()[patchi]

::
    
    kineticTheoryModels/frictionalStressModel/ReugeRahimiBoundary/ReugeRahimiBoundaryFrictionalStress.C: In member function ‘virtual Foam::tmp<Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh> > Foam::kineticTheoryModels::frictionalStressModels::ReugeRahimiBoundary::nu(const Foam::phaseModel&, const dimensionedScalar&, const dimensionedScalar&, const volScalarField&, const volSymmTensorField&, const volScalarField&, const volScalarField&) const’:
    kineticTheoryModels/frictionalStressModel/ReugeRahimiBoundary/ReugeRahimiBoundaryFrictionalStress.C:182:52: error: invalid types ‘<unresolved overloaded function type>[Foam::label {aka int}]’ for array subscript
                        sqrt( (D.boundaryField[patchi] && D.boundaryField[patchi]) + (Theta.boundaryField[patchi]/sqr(da.boundaryField[patchi]))  )
                                                        ^
    kineticTheoryModels/frictionalStressModel/ReugeRahimiBoundary/ReugeRahimiBoundaryFrictionalStress.C:182:79: error: invalid types ‘<unresolved overloaded function type>[Foam::label {aka int}]’ for array subscript
                        sqrt( (D.boundaryField[patchi] && D.boundaryField[patchi]) + (Theta.boundaryField[patchi]/sqr(da.boundaryField[patchi]))  )
                                                                                ^
    kineticTheoryModels/frictionalStressModel/ReugeRahimiBoundary/ReugeRahimiBoundaryFrictionalStress.C:182:111: error: invalid types ‘<unresolved overloaded function type>[Foam::label {aka int}]’ for array subscript
                        sqrt( (D.boundaryField[patchi] && D.boundaryField[patchi]) + (Theta.boundaryField[patchi]/sqr(da.boundaryField[patchi]))  )
                                                                                                                ^
    kineticTheoryModels/frictionalStressModel/ReugeRahimiBoundary/ReugeRahimiBoundaryFrictionalStress.C:182:140: error: invalid types ‘<unresolved overloaded function type>[Foam::label {aka int}]’ for array subscript
                        sqrt( (D.boundaryField[patchi] && D.boundaryField[patchi]) + (Theta.boundaryField[patchi]/sqr(da.boundaryField[patchi]))  )
                                                                                                                                                ^
    kineticTheoryModels/frictionalStressModel/ReugeRahimiBoundary/ReugeRahimiBoundaryFrictionalStress.C:170:27: warning: unused variable ‘U’ [-Wunused-variable]
        const volVectorField& U = phase.U();
    
Error 23 - incomplete type
--------------------------    

#include "twoPhaseSystem.H"
because  refCast<const twoPhaseSystem>(phase.fluid()).drag(phase).K() was needed

::
    
    kineticTheoryModels/viscosityModel/Srivastava/SrivastavaViscosity.C:81:53: error: invalid use of incomplete type ‘const class Foam::twoPhaseSystem’
             refCast<const twoPhaseSystem>(phase.fluid()).drag(phase).K()

Error 24 - not declared
----------------------- 

#include "twoPhaseSystem.H"
because  refCast<const twoPhaseSystem>(phase.fluid()).drag(phase).K() was needed

::

    error: ‘phase’ was not declared in this scope
         refCast<const twoPhaseSystem>(phase.fluid()).drag(phase).K()


In file included from kineticTheoryModels/conductivityModel/conductivityModel/conductivityModel.C:26:0:
kineticTheoryModels/conductivityModel/conductivityModel/conductivityModel.H:115:19: error: ‘phaseModel’ does not name a type
             const phaseModel& phase


Error 25 - file exists
----------------------

ln: failed to create symbolic link './SrivastavaCriticalStateHypothesisFrictionalStress.C': File exists
ln: failed to create symbolic link './SrivastavaCriticalStateHypothesisFrictionalStress.H': File exists
ln: failed to create symbolic link './newFrictionalStressModel.C': File exists
ln: failed to create symbolic link './frictionalStressModel.H': File exists
ln: failed to create symbolic link './frictionalStressModel.C': File exists


Error 26 - not a member
-----------------------

In file included from phaseCompressibleTurbulenceModels.C:109:0:
lnInclude/kineticTheoryModel.H:102:21: error: ‘frictionalStressModel’ is not a member of ‘Foam::kineticTheoryModels’
             autoPtr<kineticTheoryModels::frictionalStressModel>
                     ^
lnInclude/kineticTheoryModel.H:102:21: note: suggested alternative:
In file included from lnInclude/kineticTheoryModel.H:58:0,
                 from phaseCompressibleTurbulenceModels.C:109:
lnInclude/frictionalStressModel.H:52:7: note:   ‘Foam::frictionModel::frictionalStressModel’
 class frictionalStressModel
       ^
In file included from phaseCompressibleTurbulenceModels.C:109:0:
lnInclude/kineticTheoryModel.H:102:21: error: ‘frictionalStressModel’ is not a member of ‘Foam::kineticTheoryModels’
             autoPtr<kineticTheoryModels::frictionalStressModel>
                     ^
lnInclude/kineticTheoryModel.H:102:21: note: suggested alternative:
In file included from lnInclude/kineticTheoryModel.H:58:0,
                 from phaseCompressibleTurbulenceModels.C:109:
lnInclude/frictionalStressModel.H:52:7: note:   ‘Foam::frictionModel::frictionalStressModel’
 class frictionalStressModel
       ^
In file included from phaseCompressibleTurbulenceModels.C:109:0:
lnInclude/kineticTheoryModel.H:102:63: error: template argument 1 is invalid
             autoPtr<kineticTheoryModels::frictionalStressModel>
                                                               ^
/opt/openfoam5/wmake/rules/General/transform:25: recipe for target 'Make/linux64GccDPInt32Opt/phaseCompressibleTurbulenceModels.o' failed
make: *** [Make/linux64GccDPInt32Opt/phaseCompressibleTurbulenceModels.o] Error 1

Error 27 - undefined reference
------------------------------

/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsFriction.so: undefined reference to `Foam::frictionModels::frictionalStressTensorModel::dictionaryConstructorTablePtr_'
/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsFriction.so: undefined reference to `Foam::frictionModel::typeName'
/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsFriction.so: undefined reference to `Foam::frictionModels::frictionalStressTensorModel::~frictionalStressTensorModel()'
/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsFriction.so: undefined reference to `Foam::frictionModels::frictionalStressTensorModel::New(Foam::dictionary const&)'
/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsFriction.so: undefined reference to `Foam::frictionModels::frictionalStressTensorModel::frictionalStressTensorModel(Foam::dictionary const&)'
/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsFriction.so: undefined reference to `Foam::frictionModels::frictionalStressTensorModel::constructdictionaryConstructorTables()'
/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsFriction.so: undefined reference to `typeinfo for Foam::frictionModels::frictionalStressTensorModel'
/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/lib/libphaseCompressibleTurbulenceModelsFriction.so: undefined reference to `Foam::frictionModels::frictionalStressTensorModel::destroydictionaryConstructorTables()'
collect2: error: ld returned 1 exit status
/opt/openfoam5/wmake/makefiles/general:140: recipe for target '/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/bin/twoPhaseEulerFrictionFoam' failed
make: *** [/home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/bin/twoPhaseEulerFrictionFoam] Error 1

move frictionModels outside phaseCompressibleTurbulenceModels

Add this to twoPhaseEulerFrictionFoam/Make/options

::

    EXE_INC = \
    -I$(WM_PROJECT_USER_DIR)/applications/solvers/multiphase/twoPhaseEulerFrictionFoam/frictionModels/lnInclude \

    EXE_LIBS = \
    -lcompesssibleEulerianFrictionModelsFriction \

Add this to twoPhaseEulerFrictionFoam/frictionModels/Make/files

::

    frictionModels/frictionalStressTensorModel/frictionalStressTensorModel.C
    frictionModels/frictionalStressTensorModel/newfrictionalStressTensorModel.C

    frictionModels/SrivastavaSundaresan/SrivastavaSundaresan.C


    LIB = $(FOAM_USER_LIBBIN)/libcompressibleEulerianFrictionModelsFriction

Add this to twoPhaseEulerFrictionFoam/frictionModels/Make/options

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

Error 28 - frictionModel not added
----------------------------------

-lm -o /home/apr207/OpenFOAM/apr207-5.0/platforms/linux64GccDPInt32Opt/bin/twoPhaseEulerFrictionFoam
/usr/bin/ld: cannot find -lcompesssibleEulerianFrictionModelsFriction

add this to Allwmake

::

    wmake $targetType frictionModels

Or you missspelled the library!!!


add this to Allwclean

::
    
    wclean libso frictionModels
    
    
    
Error 29 - tolerance in mesh was too large
------------------------------------------

::

    Checking geometry...
        Overall domain bounding box (-0.02537582483 0 -0.001107932439) (0 0.085 0.001107932439)
        Mesh has 2 geometric (non-empty/wedge) directions (1 1 0)
        Mesh has 3 solution (non-empty) directions (1 1 1)
        Wedge symmetry_one with angle 2.512938305 degrees
    ***Wedge patch symmetry_one not planar. Point (-0.02040118265 0.07227942308 0.0008907348711) is not in patch plane by 3.001854803e-08 metre.
        Boundary openness (-4.16314945e-15 -4.881437822e-18 -1.538499505e-13) OK.
        Max cell openness = 2.845695652e-16 OK.
        Max aspect ratio = 4.406563772 OK.
        Minimum face area = 2.723616961e-09. Maximum face area = 5.530524205e-07.  Face area magnitudes OK.
        Min volume = 6.704287904e-13. Max volume = 1.367974692e-10.  Total volume = 8.947877738e-07.  Cell volumes OK.
        Mesh non-orthogonality Max: 66.93583347 average: 30.65688491
        Non-orthogonality check OK.
        Face pyramids OK.
        Max skewness = 0.8473970306 OK.
        Coupled point location match (average 0) OK.

    Failed 1 mesh checks.

https://bugs.openfoam.org/view.php?id=2126    
    
Change all tolerances in pointwise to 1e-16 (as SMALL = 1e-15 in OpenFOAM)    

Model size = 0.085
Node = 1e-16
Connector = 1e-16
Grid point = 1e-16


Error 30 - spaces in file path
------------------------------

fileName::stripInvalid() called for invalid fileName /media/apr207/MyBookDuo/0_Exeter_Desktop/home/OpenFOAM/apr207-5.0/1001_srivastava_corrected
    For debug level (= 2) > 1 this is considered fatal
    
Install GParted

Select /dev/sdc/ My Book Duo

Right click and label it My_Book_Duo

Click the tick to apply all operations



Error 31 - Cant do minimum of a matrix
--------------------------------------

Foam::tmp<Foam::fvMatrix<Foam::Vector<double> > > doesnt contain member min

::

    fvVectorMatrix dDRR(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

    dDRR =
    (
        phase2.turbulence().divDevRhoReff(U2)
    );

    Info<< "min(diagonal divDevRhoReff): " << min(dDRR.diag()) << endl;
    
    Info<< "min(upper divDevRhoReff): " << min(dDRR.upper()) << endl;
    Info<< "min(lower divDevRhoReff): " << min(dDRR.lower()) << endl;
    
    
    
Error 32 - No rule
------------------


filename or directory is wrong

::

    make: *** No rule to make target 'Make/linux64GccDPInt32Opt/kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleVelocity/JohnsonJacksonParticleVelocityFvPatchVectorField.C.dep', needed by 'Make/linux64GccDPInt32Opt/kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleVelocity/JohnsonJacksonParticleVelocityFvPatchVectorField.o'. Stop.
    

Error 33 - volVectorField is not the base of const Foam::phaseModel
-------------------------------------------------------------------
    
::
    
    kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleVelocity/JohnsonJacksonParticleVelocityFvPatchVectorField.C: In member function ‘virtual void Foam::JohnsonJacksonParticleVelocityFvPatchVectorField::updateCoeffs()’:
    kineticTheoryModels/derivedFvPatchFields/JohnsonJacksonParticleVelocity/JohnsonJacksonParticleVelocityFvPatchVectorField.C:194:36: error: ‘Foam::volVectorField {aka Foam::GeometricField<Foam::Vector<double>, Foam::fvPatchField, Foam::volMesh>}’ is not a base of ‘const Foam::phaseModel’
                phased.volVectorField::name()   
    
    
    
Error 34 - courant number is too high
-------------------------------------

Use adjustable timestep


Error 35 - cant find controlDict
--------------------------------

::
    --> FOAM FATAL ERROR:
    cannot find file "/scratch/apr207/01_06_2018_28_processors/processor0/system/controlDict"

Error 36 - output files are at the wrong times
----------------------------------------------

::

    writeControl adjustableTimeStep;
    
Error 37 - Floating point error with no exact cause
---------------------------------------------------    
    
::    
    
    #0 Foam::error::PrintStack(Foam::Ostream&) in "/opt/openfoam171/lib/linux64GccDPOpt/libOpenFOAM.so"
    #1 Foam::sigFpe::sigFpeHandler(int) in "/opt/openfoam171/lib/linux64GccDPOpt/libOpenFOAM.so"
    #2 in "/lib/libc.so.6"
    #3 Foam::GAMGSolver::scalingFactor(Foam::Field<double >&, Foam::Field<double> const&, Foam::Field<double> const&, Foam::Field<double> const&) const in "/opt/openfoam171/lib/linux64GccDPOpt/libOpenFOAM.so"
    #4 Foam::GAMGSolver::scalingFactor(Foam::Field<double >&, Foam::lduMatrix const&, Foam::Field<double>&, Foam::FieldField<Foam::Field, double> const&, Foam::UPtrList<Foam::lduInterfaceField const> const&, Foam::Field<double> const&, unsigned char) const in "/opt/openfoam171/lib/linux64GccDPOpt/libOpenFOAM.so"
    #5 Foam::GAMGSolver::Vcycle(Foam::PtrList<Foam::lduMa trix::smoother> const&, Foam::Field<double>&, Foam::Field<double> const&, Foam::Field<double>&, Foam::Field<double>&, Foam::Field<double>&, Foam::PtrList<Foam::Field<double> >&, Foam::PtrList<Foam::Field<double> >&, unsigned char) const in "/opt/openfoam171/lib/linux64GccDPOpt/libOpenFOAM.so"
    #6 Foam::GAMGSolver::solve(Foam::Field<double>&, Foam::Field<double> const&, unsigned char) const in "/opt/openfoam171/lib/linux64GccDPOpt/libOpenFOAM.so"
    #7 Foam::fvMatrix<double>::solve(Foam::dictionary const&) in "/opt/openfoam171/lib/linux64GccDPOpt/libfiniteVolume.so"
    #8 
    in "/opt/openfoam171/applications/bin/linux64GccDPOpt/buoyantSimpleFoam"
    #9 __libc_start_main in "/lib/libc.so.6"
    #10     
    
You must re-number mesh    
    
::

    renumberMesh
    
 Error 38 - Read only file system
 --------------------------------
  
 ::
    
    /usr/bin/ld: cannot open output file /usr/local/OpenFOAM/OpenFOAM-2.2.2/platforms/linux64GccDPOpt/bin/twoPhaseEulerGravityVibrationFoam: Read-only file system  
    

change FOAM_APPBIN into FOAM_USER_APPBIN in the Make/files file
    
    
Error 39 - MULES limiter
------------------------    
    
remove MULES    
    
::    
    
    Time = 0.1122
    [0] #0  Foam::error::printStack(Foam::Ostream&) at ??:?
    [0] #1  Foam::sigFpe::sigHandler(int) at ??:?
    [0] #2   in "/lib64/libc.so.6"
    [0] #3  void Foam::MULES::limiter<Foam::geometricOneField, Foam::zeroField, Foam::zeroField>(Foam::Field<double>&, Foam::geometricOneField const&, Foa$
    [0] #4  void Foam::MULES::limit<Foam::geometricOneField, Foam::zeroField, Foam::zeroField>(Foam::geometricOneField const&, Foam::GeometricField<double$
    [0] #5  Foam::MULES::explicitSolve(Foam::GeometricField<double, Foam::fvPatchField, Foam::volMesh>&, Foam::GeometricField<double, Foam::fvsPatchField,$
    [0] #6
    [0]  at ??:?
    [0] #7  __libc_start_main in "/lib64/libc.so.6"
    [0] #8
    [0]  at ??:?
    [emps-rodrigo:71579] *** Process received signal ***
    [emps-rodrigo:71579] Signal: Floating point exception (8)
    
    
::

    //
    // APR: removed MULES
    //

    //
    // APR: added 2.1.0 alphaEqn
    //

    fvScalarMatrix alpha1Eqn
    (
         fvm::ddt(alpha1)
       + fvm::div(phic, alpha1, alphaScheme)
       + fvm::div(-fvc::flux(-phir, alpha2, alpharScheme), alpha1, alpharScheme)
    );

    if (g0.value() > 0.0)
    {

        surfaceScalarField alpha1f(fvc::interpolate(alpha1));

        ppMagf =
            rAU1f/(alpha1f + scalar(0.0001))
           *(g0/rho1)*min(exp(preAlphaExp*(alpha1f - alphaMax)), expMax);

        alpha1Eqn -= fvm::laplacian
        (
            (fvc::interpolate(alpha1) + scalar(0.0001))*ppMagf,
            alpha1,
            "laplacian(alpha1PpMag,alpha1)"
        );
    }

    //
    // APR: added 2.1.0 alphaEqn
    //

    
    
Error 40 - exp floating point
-----------------------------    
 
interpolate whole of ppMagf
 
    
::    

    Time = 0.0914
    [0] #0  Foam::error::printStack(Foam::Ostream&)[7] #0  Foam::error::printStack(Foam::Ostream&)[6] #0  Foam::error::printStack(Foam::Ostream&)[8] #0  F$
    [7] #1  Foam::sigFpe::sigHandler(int) at ??:?
    [3] #2   at ??:?
    [5] #2   in "/lib64/libc.so.6"
    [0] #4  exp in "/lib64/libm.so.6"
    [0] #5  Foam::exp(Foam::Field<double>&, Foam::UList<double> const&) in "/lib64/libm.so.6"
    [emps-rodrigo:88525] *** Process received signal ***
    [emps-rodrigo:88525] Signal: Floating point exception (8)
    [emps-rodrigo:88525] Signal code:  (-6)
    [emps-rodrigo:88525] Failing at address: 0x76feb1000159cd    
  
    
::
    
    //surfaceScalarField alpha1f(fvc::interpolate(alpha1));

    //
    // APR: interpolate all for ppMagf
    //
        ppMagf = rAU1f*fvc::interpolate
        (
            (1.0/(rho1*(alpha1 + scalar(0.0001))))
           *g0*min(exp(preAlphaExp*(alpha1 - alphaMax)), expMax)
        );
    //
    // APR: interpolate all for ppMagf
    //
    
    
    
    
Error 41 - exp floating point
----------------------------- 
    
::

    Time = 0.2461
    [2] #0  Foam::error::printStack(Foam::Ostream&)[0] #0  Foam::error::printStack(Foam::Ostream&)[1] #0  Foam::error::printStack(Foam::Ostream&) at ??:?
    [2] #1  Foam::sigFpe::sigHandler(int) at ??:?
    [2] #2   at ??:?
    [1] #2   in "/lib/x86_64-linux-gnu/libc.so.6"
    [2] #4  exp in "/lib/x86_64-linux-gnu/libm.so.6"
    [2] #5  Foam::exp(Foam::Field<double>&, Foam::UList<double> const&) in "/lib/x86_64-linux-gnu/libm.so.6"
    [emps-andrew:12855] *** Process received signal ***
    [emps-andrew:12855] Signal: Floating point exception (8)
    [emps-andrew:12855] Signal code:  (-6)
    [emps-andrew:12855] Failing at address: 0x3e800003237   
    
    
    
Error 42 - Simulation just stops with axi-symmetric wedge BCs
-------------------------------------------------------------

Recommend wedge angle is less than 5 degrees, i.e. 1 degree (not the problem)

compare working state log.twoPhaseEulerSrivastavaFoam (non-wedge) with non-working state (wedge).

Found that it was getting stuck on the line set, so I removed #include singleGraph


Error 43 - High air velocities at boundary
------------------------------------------

Air is escaping from the domain and cannot escape, because outlet air meets a wall

So replace outletInlet with inletOutlet at the inlet


Error 44 - Cannot find p file
-----------------------------

::

    [0] --> FOAM FATAL ERROR:
    [0] cannot find file "/scratch/apr207/23_inlet_bc_and_momentum_predictor/processor0/0/p"
    [0]
    [0]     From function virtual Foam::autoPtr<Foam::ISstream> Foam::fileOperations::uncollatedFileOperation::readStream($
    [0]     in file global/fileOperations/uncollatedFileOperation/uncollatedFileOperation.C at line 505.

Look in log.decomposePar:

::

    --> FOAM FATAL IO ERROR:
    keyword inletValue is undefined in dictionary "/scratch/apr207/23_inlet_bc_and_momentum_predictor/0/U.air.boundaryField
    
