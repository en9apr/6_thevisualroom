====================
Installation History
====================

.. contents::
   :local:

Remove message
==============

::

    touch ~/.sudo_as_admin_successful

OpenFOAM 1612+
==============
http://openfoam.com/download/install-source.php

http://openfoam.com/code/build-guide.php

OpenFOAM 4.1
============

https://openfoam.org/download/4-1-ubuntu/

OpenFOAM 5.0
============

https://openfoam.org/download/5-0-ubuntu/

OpenFOAM Extend 4.0
===================

https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.0

OpenFOAM 2.4.0
==============

https://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-2.4.0/Ubuntu

How to create extension .foam with paraFoam -touch?
===================================================

::

    $ nano OpenFOAM/OpenFOAM-2.4.0/bin/paraFoam

    # reader extension
    extension=foam


Kate
====

::

    sudo apt-get update
    sudo apt-get install kate

    Show line numbers:

    Settings -> configure kate -> appearance -> borders 

Dropbox
=======

Ubuntu Software Centre then Daemon

::

    sudo apt-get update  //update package list - need to log out and log in to restart Nautilius

Texmaker
========

::

    sudo apt-get update
    sudo apt-get install texlive-full
    sudo apt-get install texmaker
    http://www.artfiles.org/openoffice.org/contrib/dictionaries/en_GB.zip
    Unzip the archive
    cd /usr/share/hunspell
    sudo cp -rv /home/apr207/Downloads/en_GB.dic .
    sudo cp -rv /home/apr207/Downloads/en_GB.aff .

SSH
===

::

    sudo apt-get install openssh-server
    from laptop: ssh -XC user@IP address

    Start ssh: sudo service ssh start
    Stop ssh: service ssh stop

Gimp
====

Ubuntu Software Centre

(Evince gave an error)

Pointwise
=========

::

    Downloaded pw-V18.0R2-linux_x86_64-jre.sh
    sh pw-V18.0R2-linux_x86_64-jre.sh
    Next
    Accept
    /home/apr207/Pointwise/PointwiseV18.0R2 is installation directory
    Run Pointwise
    (Won't find license)
    Specify license server
    server: emps-pointwise
    port: 2385

    Added this to .bashrc:
    # Pointwise
    alias pointwise="/home/apr207/Pointwise/PointwiseV18.0R2/pointwise"

Tree
====

::

    sudo apt-get install tree


LIGGGHTS
========

::

    ## [optional] 1. Install Voro++ 0.4.x by compiling
    sudo apt install g++
    cd ~
    wget http://math.lbl.gov/voro++/download/dir/voro++-0.4.6.tar.gz
    tar -zxvf voro++-0.4*.tar.gz
    cd $HOME/voro++-0.4.*
    make all
    sudo make install
    which voro++
    # /usr/local/bin/voro++ <- comes up
    #
    ## 2 Install Liggghts 3.x by compiling
    ## 2.1 Install Packages
    sudo apt install git libvtk5-dev libeigen2-dev openmpi-bin openmpi-doc libopenmpi-dev
    which mpirun
    #/usr/bin/mpirun comes up
    which mpic++
    #/usr/bin/mpic++ comes up
    ## 2.2 Get Liggghts via Git
    cd ~
    git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC $HOME/LIGGGHTS-PUBLIC3.6.0
    ## 2.3 Compiling Liggghts with VORONOI and jpg, png support
    cd $HOME/LIGGGHTS-PUBLIC3.6.0/src
    # [optional] if you want to use voro++ in LIGGGHTS: 
    make yes-voronoi
    # [optional] if you need extra packages install with sudo make yes-packagename
    gedit /$HOME/LIGGGHTS-PUBLIC3.6.0/src/MAKE/Makefile.ubuntuVTK
    add "-DLAMMPS_JPEG -DLAMMPS_PNG" in line 32
    add "-I/usr/include" in line 63
    add "-ljpeg -lpng" in line 65
    change line 73 to "VTK_INC = -I/usr/include/vtk-5.10"
    change line 74 to "VTK_LIB = -lvtkCommon -lvtkFiltering -lvtkIO -lvtkParallel -lvtkGraphics"
    # save and close gedit
    # [optional] changes on the source code
    make clean-all
    make ubuntuVTK
    # create system wide shortcut liggghts360 for compiled binary (I used to have different versions parallel)
    sudo ln -s /$HOME/LIGGGHTS-PUBLIC3.6.0/src/lmp_ubuntuVTK /usr/bin/liggghts360
    liggghts360

    # Liggghts comes up, telling version, compiling date etc., press Ctr+d to quit
    #
    ## [optional] 3. Install LPP for post processing (converts LIGGGHTS Dumps to vtk-files) - doesn't work?
    cd ~
    sudo apt-get install python-numpy
    #already newest version
    sudo git clone https://github.com/CFDEMproject/LPP.git $HOME/LPP
    ./install
    gedit ~/.bashrc
    #add: 
    alias lpp="python $HOME/LPP/src/lpp.py"
    #sudo chown -R andrew:andrew LPP
    #open new Terminal: lpp
    #
    ## [optional] 4. Install Syntax Highlighting for xed (gedit)
    cd ~
    wget https://www.dropbox.com/s/78elqj4i2dn52wt/liggghts3.lang
    sudo mv liggghts3.lang /usr/share/gtksourceview-3.0/language-specs
    #
    ## [optional] 5. Install ParaView 5.0.1
    sudo apt-get install paraview
    # already newest version
    ## [optional] 6. Install GNUplot 5.0.3
    sudo apt-get install gnuplot-x11
    gnuplot
    plot sin(x)
    # window with sin graph comes up, press Ctr+d to quit.
    #
    ## [optional] 7. Install Povray 3.7.1
    # Alternative: http://www.conoce3000.com/html/espaniol/Apuntes/2014-06-20-CompilarInstalarPOV-Ray37LinuxMintCinnamon64bitsCompilarInstalarPOV-Ray37LinuxMintCinnamon64bits.php?Arch=20
    cd ~
    sudo apt-get install autoconf automake libboost-all-dev libboost-dev libopenexr-dev libsdl-dev zlib1g-dev libpng-dev libjpeg-dev libtiff-dev
    git clone https://github.com/POV-Ray/povray.git $HOME/POV-Ray3.7
    cd $HOME/POV-Ray3.7/unix
    ./prebuild.sh
    cd ..
    ./configure COMPILED_BY="andrew"
    make check
    # windows with cup and cookies comes up, click picture
    sudo make install
    # done!

    # Test LIGGHTS

    Copy all tutorials to LIGGGHTS_User folder

    cd $HOME/LIGGGHTS_User/chute_wear

    liggghts360 < in.chute_wear

    cd post

    lpp dump*.chute

VIM
===

::

    sudo apt-get install vim

INKSCAPE
========

::

    sudo apt-get install inkscape

LAMMPS
======

::

    sudo add-apt-repository ppa:gladky-anton/lammps
    sudo apt-get update 

    sudo apt-get install lammps-daily 

    lammps-daily < in.lj 

    sudo apt-get update 

    sudo apt-get install lammps-daily-doc 


CFMesh
======

::

    Downloaded CFMesh: 

    https://sourceforge.net/projects/cfmesh/

    Copied cfMesh-v1.1.2 to /home

    Set environment to $ two (OpenFOAM version 2.4 - as this is installed on Rodrigo and is possible on Callisto - no instruction for 2.3.1)

    $ ./Allwmake

    Copy the tutorial files:

    $ cp -rf tutorials $FOAM_RUN


    $ file file.stl

    file.stl: ASCII text

Mesh Lab
========

Ubuntu Software Centre

FreeCAD
=======

Ubuntu Software Centre


Pizza.py
========

::

    Downloaded Pizza.py from 

    https://sourceforge.net/projects/pizza-py/?source=typ_redirect

    Extracted the file
    Copied to /home

    Added this to .bashrc

    #Pizza.py
    alias pizza="python -i $HOME/pizza/src/pizza.py"

Installed Redshift
==================

Ubuntu Software Centre



Install Heekes CAD
==================

::

    Add heekscnc-devel PPA to your repositories list:

    ### sudo add-apt-repository ppa:neomilium/heekscnc-devel ### maybe not the development version

    sudo add-apt-repository ppa:neomilium/cam
    
    Update your packages list:

    sudo apt-get update

Install HeeksCNC
================

::

    sudo apt-get install heekscnc

Install Mendeley
================

Ubuntu Software Centre
64 bit version of Mendeley: https://www.mendeley.com/download-mendeley-desktop/ubuntu/instructions/
Install automatically

Install Pinta
=============

Ubuntu Software Centre

Install Curl
============

::

    sudo apt install curl

Check for Updates
=================

To check for updates: Start > Software Updater

Check packages have been updated
================================

::

    /usr/lib/update-notifier/apt-check -p

Install Scipy stack
===================

::

    $ sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose

Remove gedit
============

Ubuntu Software Centre

Terminator
==========

Ubuntu Software Centre

Sublime
=======

::

    sudo rpm -v --import https://download.sublimetext.com/sublimehq-rpm-pub.gpg

    sudo dnf config-manager --add-repo https://download.sublimetext.com/rpm/stable/x86_64/sublime-text.repo

    sudo dnf install sublime-text

    Edit sublime_text.desktop

    [Desktop Entry]
    Encoding=UTF-8
    Version=1.0
    Type=Application
    Name=Sublime Text
    Icon=sublime_text.png
    Path=/
    Exec=/opt/sublime_text/sublime_text %f
    StartupNotify=false
    StartupWMClass=Sublime_text
    OnlyShowIn=Unity;
    X-UnityGenerated=true


OpenFOAM 2.1.x
==============

::

    FROM: http://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-2.1.x/Ubuntu

    Up to step 7.

    FROM: https://www.cfd-online.com/Forums/openfoam-installation/168746-problems-installing-openfoam-2-4-0-ubuntu-16-04-a.html


    #Go into OpenFOAM's main source folder
    cd $WM_PROJECT_DIR
     
    #Change how the flex version is checked
    find src applications -name "*.L" -type f | xargs sed -i -e 's=\(YY\_FLEX\_SUBMINOR\_VERSION\)=YY_FLEX_MINOR_VERSION < 6 \&\& \1='

    #Still better be certain that the correct Qt version is being used
    export QT_SELECT=qt4

    #Back to step 8

    cd $WM_PROJECT_DIR

    # This next command will take a while... somewhere between 30 minutes to 3-6 hours.
    ./Allwmake > make.log 2>&1
     
    #Run it a second time for getting a summary of the installation
    ./Allwmake > make.log 2>&1


OpenFOAM 2.1.0
==============

::

    Till step 4:

    http://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-2.1.x/Ubuntu#Ubuntu_14.04

    sudo -s

    apt-get update

    apt-get install build-essential cmake flex bison zlib1g-dev qt4-dev-tools libqt4-dev gnuplot libreadline-dev \
    libncurses-dev libxt-dev libopenmpi-dev openmpi-bin git-core gcc-4.7 g++-4.7

    apt-get install libscotch-dev

    exit

    https://openfoam.org/download/2-1-0-source/

    download OpenFOAM-2.1.0 and ThirdParty-2.1.0

    cd ThirdParty-2.1.0
    mkdir download
    wget -P download http://www.paraview.org/files/v3.12/ParaView-3.12.0.tar.gz
    wget -P download https://gforge.inria.fr/frs/download.php/28043/scotch_5.1.11.tar.gz
    tar -xzf download/ParaView-3.12.0.tar.gz
    tar -xzf download/scotch_5.1.11.tar.gz
    cd ..

    Now step 6:

    uname -m
    sed -i -e 's/gcc/\$(WM_CC)/' OpenFOAM-2.1.0/wmake/rules/linux64Gcc/c
    sed -i -e 's/g++/\$(WM_CXX)/' OpenFOAM-2.1.0/wmake/rules/linux64Gcc/c++
    source $HOME/OpenFOAM/OpenFOAM-2.1.0/etc/bashrc WM_NCOMPPROCS=8 WM_MPLIB=SYSTEMOPENMPI
    export WM_CC='gcc-4.7'
    export WM_CXX='g++-4.7'

    FULL_SETTINGS="$FOAM_SETTINGS; export WM_CC=gcc-4.7; export WM_CXX=g++-4.7"
    echo "alias of210='source \$HOME/OpenFOAM/OpenFOAM-2.1.0/etc/bashrc $FULL_SETTINGS'" >> $HOME/.bashrc
    unset FULL_SETTINGS

    cd OpenFOAM-2.1.0
    ./Allwmake > log.make 2>&1
    ./Allwmake > log.make 2>&1

    FROM: http://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-2.1.x/Ubuntu

    Up to step 7.

    FROM: https://www.cfd-online.com/Forums/openfoam-installation/168746-problems-installing-openfoam-2-4-0-ubuntu-16-04-a.html

    #Go into OpenFOAM's main source folder
    cd $WM_PROJECT_DIR
     
    #Change how the flex version is checked
    find src applications -name "*.L" -type f | xargs sed -i -e 's=\(YY\_FLEX\_SUBMINOR\_VERSION\)=YY_FLEX_MINOR_VERSION < 6 \&\& \1='

    #Still better be certain that the correct Qt version is being used
    export QT_SELECT=qt4

    #Back to step 8

    cd $WM_PROJECT_DIR

    # This next command will take a while... somewhere between 30 minutes to 3-6 hours.
    ./Allwmake > make.log 2>&1
     
    #Run it a second time for getting a summary of the installation
    ./Allwmake > make.log 2>&1


Generate OpenFOAM 2.2.2 doxygen
===============================

::

    sudo apt-get install doxygen graphviz

    ./Allmake

OpenFOAM 2.3.0
==============

::

    # The following lines shouldn't do anything:

    sudo -s

    apt-get update

    apt-get install build-essential cmake flex bison zlib1g-dev qt4-dev-tools libqt4-dev libqtwebkit-dev gnuplot \
    libreadline-dev libncurses5-dev libxt-dev libopenmpi-dev openmpi-bin libboost-system-dev libboost-thread-dev libgmp-dev \
    libmpfr-dev python python-dev

    apt-get install libglu1-mesa-dev libqt4-opengl-dev

    exit

    # download OpenFOAM-2.3.0 and ThirdParty-2.3.0:

    cd ~

    cd OpenFOAM

    wget "http://downloads.sourceforge.net/foam/OpenFOAM-2.3.0.tgz?use_mirror=mesh" -O OpenFOAM-2.3.0.tgz

    wget "http://downloads.sourceforge.net/foam/ThirdParty-2.3.0.tgz?use_mirror=mesh" -O ThirdParty-2.3.0.tgz
     
    tar -xzf OpenFOAM-2.3.0.tgz 
    tar -xzf ThirdParty-2.3.0.tgz

    # Symbolic links:

    ln -s /usr/bin/mpicc.openmpi OpenFOAM-2.3.0/bin/mpicc
    ln -s /usr/bin/mpirun.openmpi OpenFOAM-2.3.0/bin/mpirun

    uname -m

    source $HOME/OpenFOAM/OpenFOAM-2.3.0/etc/bashrc WM_NCOMPPROCS=8 WM_MPLIB=SYSTEMOPENMPI

    echo "alias of230='source \$HOME/OpenFOAM/OpenFOAM-2.3.0/etc/bashrc $FOAM_SETTINGS'" >> $HOME/.bashrc

    # Open new terminal:

    of230

    # Paraview

    cd $WM_THIRD_PARTY_DIR
     
    #make very certain that the correct Qt version is being used, by running this command:
    export QT_SELECT=qt4
     
    # This next command will take a while... somewhere between 5 minutes to 30 minutes.
    ./Allwmake > log.make 2>&1
     
    #update the shell environment
    wmSET $FOAM_SETTINGS

    #Go into OpenFOAM's main source folder
    cd $WM_PROJECT_DIR
     
    #Change how the flex version is checked
    find src applications -name "*.L" -type f | xargs sed -i -e 's=\(YY\_FLEX\_SUBMINOR\_VERSION\)=YY_FLEX_MINOR_VERSION < 6 \&\& \1='

    #Still better be certain that the correct Qt version is being used
    export QT_SELECT=qt4

    #Go into OpenFOAM's main source folder
    cd $WM_PROJECT_DIR
     
    #Still better be certain that the correct Qt version is being used
    export QT_SELECT=qt4
     
    # This next command will take a while... somewhere between 30 minutes to 3-6 hours.
    ./Allwmake > log.make 2>&1
     
    #Run it a second time for getting a summary of the installation
    ./Allwmake > log.make 2>&1

    icoFoam -help


OpenFOAM 2.2.x
==============

::

    http://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-3.0.1/Ubuntu#Ubuntu_16.04

    sudo -s

    apt-get update

    apt-get install build-essential flex bison zlib1g-dev qt4-dev-tools libqt4-dev gnuplot libreadline-dev \
    libncurses-dev libxt-dev libopenmpi-dev openmpi-bin

    exit

    cd OpenFOAM

    download from github OpenFOAM-2.2.x-master and ThirdParty-2.2.x-master

    unzip using archive manager

    remove master from the folder names

    uname -m

    source $HOME/OpenFOAM/OpenFOAM-2.2.x/etc/bashrc WM_NCOMPPROCS=2 WM_MPLIB=SYSTEMOPENMPI

    echo "alias of22x='source \$HOME/OpenFOAM/OpenFOAM-2.2.x/etc/bashrc $FOAM_SETTINGS'" >> $HOME/.bashrc

    sed -i -e 's=-lz -lm -lrt=-Xlinker --no-as-needed -lz -lm -lrt=' \
      ThirdParty-2.2.x/etc/wmakeFiles/scotch/Makefile.inc.i686_pc_linux2.shlib-OpenFOAM-*
      
    #Go into OpenFOAM's main source folder
    cd $WM_PROJECT_DIR
     
    #Change how the flex version is checked
    find src applications -name "*.L" -type f | xargs sed -i -e 's=\(YY\_FLEX\_SUBMINOR\_VERSION\)=YY_FLEX_MINOR_VERSION < 6 \&\& \1='

    #Still better be certain that the correct Qt version is being used
    export QT_SELECT=qt4

    ./Allwmake > log.make 2>&1

    ./Allwmake > log.make 2>&1


To create a desktop shortcut for Doxygen
========================================

::

    Ubuntu software centre: install gnome panel

    $ gnome-desktop-item-edit --create-new ~/Desktop



Symbolic Link to OpenFOAM folder
================================

DON'T CREATE ~/Dropbox/OpenFOAM/

::

    $ ln -s /home/apr207/OpenFOAM/apr207-2.4.0/ ~/Dropbox/OpenFOAM/

    $ ln -s /home/apr207/OpenFOAM/apr207-2.2.2/ ~/Dropbox/OpenFOAM/

    $ ln -s /home/apr207/OpenFOAM/apr207-5.0/applications ~/Dropbox/OpenFOAM/apr207-5.0/applications
    
    If there is a mistake:

    cd ~/Dropbox/OpenFOAM/

    unlink apr207-2.4.0

DONT CREATE ~/Dropbox/Pointwise_User/

    $ ln -s /home/apr207/Pointwise_User/ ~/Dropbox/

    $ ln -s /home/apr207/LAMMPS_User/ ~/Dropbox/

    $ ln -s /home/apr207/LIGGGHTS_User/ ~/Dropbox/
    
    $ ln -s /home/apr207/blender-2.78c-linux-glibc219-x86_64/ ~/Dropbox/
    
    $ ln -s /home/apr207/cfMesh-v1.1.2/ ~/Dropbox/

Symbolic link to subfolder only
===============================

make directory apr207-5.0 on Dropbox

::
    
    $ ln -s /home/apr207/OpenFOAM/apr207-5.0/applications ~/Dropbox/OpenFOAM/apr207-5.0/
    
make directory run on dropbox    
    
    $ ln -s /home/apr207/OpenFOAM/apr207-5.0/run/twoPhaseEulerCohesionFoam ~/Dropbox/OpenFOAM/apr207-5.0/run/
    
    $ ln -s /home/apr207/OpenFOAM/apr207-5.0/run/twoPhaseEulerDenseFoam_hallflowmeter ~/Dropbox/OpenFOAM/apr207-5.0/run/
    
    $ 
    
DOS2UNIX
========

::

    sudo apt install dos2unix

Mencoder
========

Ubuntu Software Centre


Engauge Digitiser
=================

Ubuntu Software Centre


Foxit Reader
============

https://www.foxitsoftware.com/pdf-reader/

Sphinx
======

::

    sudo apt-get install python-sphinx

Filezilla
=========

Ubuntu Software Centre

Gitk
====

Ubuntu Software Centre

Ubuntu Freezes Randomly Solution
================================

https://askubuntu.com/questions/761706/ubuntu-15-10-and-16-04-keep-freezing-randomly

sudo nano /etc/default/grub
There is a line in that: GRUB_CMDLINE_LINUX_DEFAULT="quiet splash" (like this), replace with: GRUB_CMDLINE_LINUX_DEFAULT="quiet splash intel_idle.max_cstate=1"
Save it (CTRL+O)
sudo update-grub
sudo reboot

I think this makes the mouse disappear (!)


Installing OpenFOAM-5.0 on emps-rodrigo (emps-machines)
=======================================================

1) Download packages:

download from http://dl.openfoam.org/source/5-0
download from http://dl.openfoam.org/third-party/5-0
scp -r OpenFOAM-5.x-version-5.0.tar.gz apr207@emps-rodrigo:/home/links/apr207/OpenFOAM
scp -r ThirdParty-5.x-version-5.0.tar.gz apr207@emps-rodrigo:/home/links/apr207/OpenFOAM
tar -zxvf OpenFOAM-5.x-version-5.0.tar.gz
tar -zxvf ThirdParty-5.x-version-5.0.tar.gz
mv OpenFOAM-5.x-version-5.0 OpenFOAM-5.0
mv ThirdParty-5.x-version-5.0 ThirdParty-5.0

2) Don't need to install any software for compilation 

3) Set environment variables: 

cd $HOME
nano .bashrc
alias five="module load mpi; source $HOME/OpenFOAM/OpenFOAM-5.0/etc/bashrc"

4) After logging out and logging in again, install Third Party Scotch and PT Scotch

five
cd $WM_THIRD_PARTY_DIR
./Allwmake > log &

5) Install OpenFOAM

cd ../OpenFOAM-5.0/
./Allwmake > log &

However, the following error occurred during stage 5):

touch: cannot touch ‘/secamfs/userspace/staff/apr207/OpenFOAM/OpenFOAM-5.0/platforms/linux64GccDPInt32OptSYSTEMOPENMPI/src/parallel/decompose/ptscotchDecomp/using:openmpi-system’: No such file or directory
touch: cannot touch ‘/secamfs/userspace/staff/apr207/OpenFOAM/OpenFOAM-5.0/platforms/linux64GccDPInt32OptSYSTEMOPENMPI/src/parallel/decompose/ptscotchDecomp/using:scotch_6.0.3’: No such file or directory

6) Install Pstream separately

cd src/parallel
./Allwmake

7) Re-install OpenFOAM

cd ../OpenFOAM-5.0/
./Allwmake




Install traceroute
==================

See who is connected to remote server

sudo apt install traceroute

$ traceroute emps-kodaly

Install CFMESH
==============

Download from https://sourceforge.net/projects/cfmesh/files/v1.1.2/

Extract to ~/home/apr207/cfMesh-v1.1.2

It only works with version 2.4.0:

$ two
$ export WM_NCOMPPROCS=[8]
$ source /home/apr207/OpenFOAM/OpenFOAM-2.4.0/etc/bashrc
$ ./Allwmake


Installing OpenFOAM-2.1.0 on emps-rodrigo (emps-machines)
=========================================================

1) Download packages:

::

    Till step 4:

    http://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-2.1.x/Ubuntu#Ubuntu_14.04

    https://openfoam.org/download/2-1-0-source/

    download OpenFOAM-2.1.0 and ThirdParty-2.1.0

    scp -r OpenFOAM-2.1.0.tgz apr207@emps-rodrigo:/home/links/apr207/OpenFOAM
    
    scp -r ThirdParty-2.1.0.tgz apr207@emps-rodrigo:/home/links/apr207/OpenFOAM

    tar -zxvf OpenFOAM-2.1.0.tgz

    tar -zxvf ThirdParty-2.1.0.tgz


3) Set environment variables: 

::

    cd ~
    nano .bashrc
    alias twoone="module load mpi; source $HOME/OpenFOAM/OpenFOAM-2.1.0/etc/bashrc"


4) After logging out and logging in again, install Third Party Scotch and PT Scotch

::

    twoone
    cd $WM_THIRD_PARTY_DIR
    ./Allwmake > log 2>&1 &
    
5) Install OpenFOAM

::

    cd ../OpenFOAM-2.1.0/
    ./Allwmake > log 2>&1 &    
    

    
    could not open file omp.h for source file PstreamGlobals.C
    could not open file openmpi/ompi/mpi/cxx/mpicxx.h for source file PstreamGlobals.C

    
    
Vim highlighting for case files
===============================

Install pathogen and create .vimrc file

::

    mkdir -p ~/.vim/autoload ~/.vim/bundle && \
    curl -LSso ~/.vim/autoload/pathogen.vim https://tpo.pe/pathogen.vim

    vim ~/.vimrc
    
        execute pathogen#infect()
        syntax on
        filetype plugin indent on
        set t_Co=256
        
Clone the extension

::    
    
    cd $HOME/.vim/bundle  
    git clone https://bitbucket.org/shor-ty/vimextensionopenfoam.git
    
    
Check that 256 colours are present:

::

    export TERM=screen-256color   
    
    
    
Terminator color scheme (Lemur only)
====================================

Right click the terminator window and select preferences. Go to the Profile tab and add a new profile. Give it a name. Configure how you want the terminal to look. Close out of that.

When you launch terminator, launch it like so

::

    terminator --profile=profilename
    
That will load terminator with a profile with the name of profilename. Change profilename to what you called yours.


Redshift colour
===============

Edit those by adding the -t flag followed by the values (in the form of day:night)

::

    redshift-gtk -l 50.7:-3.53 -t 2500:2500

You'll have to play around with the numbers a bit to find one that works for you. The lower the number, the more red it will get - 6500 being no tint at all. To make RedShift start up when your computer does, you can easily do so by going to Startup Applications and adding a new program. Name it whatever you want and type the above command (using your settings) in the command box.

Simple Screen Recorder
======================

::

    sudo add-apt-repository ppa:maarten-baert/simplescreenrecorder
    sudo apt-get update
    sudo apt-get install simplescreenrecorder

Error - System program problem detected
=======================================

System program problem detected
Do you want to report the problem now?

See the crashes:

::

    ls -l /var/crash/
    sudo rm /var/crash/*
    sudo gedit /etc/default/apport &
    
        # set this to 0 to disable apport, or to 1 to enable it
        # you can temporarily override this with
        # sudo service apport start force_start=1
        enabled=0
    
    
How to install discord on Ubuntu
================================

Load terminator

::

    $ discord

It says to download discord
Download deb file
Open debfile
Done

::

    $ discord

link account 

Install Spyder
==============

Ubuntu software centre   
    
    
