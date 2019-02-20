=========
Swak4Foam
=========


Download and check requirements:

::

    $ two
    $ cd "$HOME/OpenFOAM/$USER-$WM_PROJECT_VERSION"
    $ git clone https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-Breeder2.0-libraries-swak4Foam.git swak4Foam
    $ cd swak4foam
    $ ln -s swakConfiguration.automatic swakConfiguration
    $ ./maintainanceScripts/compileRequirements.sh
  
In .bashrc:

::

    export PATH="/home/andrew/OpenFOAM/andrew-2.4.0/swak4Foam/privateRequirements/bin:$PATH"

Maybe open a new terminal Build swak4Foam:    
    
::

    $ two
    $ cd "$HOME/OpenFOAM/$USER-$WM_PROJECT_VERSION"
    $ wmake all

In .bashrc:

::
    
    export SWAK4FOAM_SRC="/home/andrew/OpenFOAM/andrew-2.4.0/swak4Foam/Libraries"
    
    
    
    
