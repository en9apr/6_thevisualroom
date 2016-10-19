==========================
Software and Customisation
==========================

* This is a list of software I've installed so far on the VM, plus additional customisation.

::

	# Guest Additions Installed

	ran autorun.sh on Guest Additions CD

	# OpenFOAM installed

	VERS=$(lsb_release -cs) 
	sudo sh -c "echo deb http://www.openfoam.org/download/ubuntu $VERS main > /etc/apt/sources.list.d/openfoam.list"
	sudo apt-get update
	sudo apt-get install openfoam231

	# paraview installed

	sudo apt-get install paraviewopenfoam410 

	# gedit installed and edited .bashrc to run OpenFOAM

	sudo apt-get install gedit

	# edited .bashrc

	gedit ~/.bashrc
	Added: source /opt/openfoam231/etc/bashrc
	export FOAM_RUN='/home/andrew/OpenFOAM/andrew-2.3.1/run'
	export WORK='/home/andrew/OpenFOAM/andrew-2.3.1/run/tutorials/incompressible/icoFoam/cavity'
	export DOWNLOADS='/home/andrew/Downloads'
	export SPHINX='/home/andrew/Sphinx/thevisualroom'

	# sphinx
	sudo apt-get install python-sphinx

	# removed Pictures, Documents, Videos, Music, Templates, Public

	# git
	sudo apt-get install git
	sudo apt-get install gitk               # to view git history

	# filezilla
	sudo apt-get install filezilla

	# Anaconda
	Downloaded Anaconda: docs.continuum.io/anaconda/install.html
	Change to Downloads directory
	bash Anaconda-2.2.0-Linux-x86_4.sh
	Accept licence and install to /home/andrew/anaconda

	sudo apt-get update

	# latex:
	sudo apt-get install texlive-full
	sudo apt-get install texmaker

	# kate:
	sudo apt-get install kate
	
	# emacs:
	sudo apt-get install emacs

	# tree - for viewing directory structures
	sudo apt-get install tree

	# Inkscape - better than GIMP
	sudo apt-get install inkscape

	# unzip program:
	sudo apt-get install unzip

	# Fortran related:
	sudo apt-get install gfortran #fortran
	sudo apt-get install liblapack-dev #lapack and blas
	sudo apt-get install openmpi-bin #openmpi executables
	sudo apt-get install libopenmpi-dev #openmpi libraries
	sudo apt-get install gdb #debugger for gfortran




