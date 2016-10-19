===============================
How to Create a Virtual Machine
===============================

.. contents::
   :local:

Why create a VM?
================

Advantages
----------

* Can test software installations in a sandbox before use, can revert changes if it doesn't work
* Can give the virtual machine with the software to the client (so they don't have to install anything to run the program)

Disadvantages
-------------

* Need enough RAM: need twice the amount of memory of the VMs
* Have to have enough disk space
* Have to have 64-bit CPU with virtualisation extensions

Process of Creating a VM
========================

`Process of Creating the VM for CentOS <http://www.tecmint.com/centos-7-installation/>`_

Download Disk Image
-------------------

* Download the .iso file (64 bit) e.g. xubuntu-14.04.2-desktop-amd64.iso

`Xubuntu Desktop <http://cdimages.ubuntu.com/xubuntu/releases/14.04/release/>`_

* Why Xubuntu:

 - It's very lightweight (maybe not as lightweight as I thought... as it has package manager)

* Why not CentOS?

 - It doesn't have package manager

* Download CentOS (minimal - if you know what you are doing... this will not install any software, even if you select it in the installation):

`Minimal CentOS <http://isoredirect.centos.org/centos/7/isos/x86_64/CentOS-7-x86_64-Minimal-1503-01.iso>`_

I used version 7 (although others use 6.6). The mirror I used to download it was this (the iso is called: CentOS-7-x86_64-Minimal-1503-01.iso)

`Minimal Centos Mirror <http://www.mirrorservice.org/sites/mirror.centos.org/7/isos/x86_64/CentOS-7-x86_64-Minimal-1503-01.iso>`_

* Download CentOS DVD:

`DVD Centos <http://isoredirect.centos.org/centos/7/isos/x86_64/CentOS-7-x86_64-DVD-1503-01.iso>`_

I used version 7 (although others use 6.6). The mirror I used to download it was this (the iso is called: CentOS-7-x86_64-DVD-1503-01.iso )

`DVD Centos Mirror <http://mirrors.coreix.net/centos/7/isos/x86_64/CentOS-7-x86_64-DVD-1503-01.iso>`_

Create a New VirtualBox
-----------------------

* Download and Install Oracle VirtualBox:

`Oracle VirtualBox <https://www.virtualbox.org/wiki/Downloads>`_

* Startup VirtualBox and create a new VM e.g.

 - Name: Xubuntu 14.04.2 LTS 20GB
 - Type: Linux 
 - Version: Xubuntu 64 bit

* Or:

 - CentOS-7-x86_64-DVD-1503-01
 - Type: Linux
 - Version: Red Hat (64-bit)

Create Virtual Hard Drive
-------------------------

* Memory size (half of what is installed): 6GB
* Select option: Create a virtual harddrive now
* Select option: VirtualBox Disk Image (experience is that the other options don't work any better)
* Select option: Fixed size 20GB (so that we have enough space for Xubuntu and other updates) ... Actually I think 50GB is better for having space to do things.
* Hit "Create" (might take about 10 minutes or so... will take 30min if it's 50GB)

You now have a virtual harddrive

Install Operating System - Ubuntu
---------------------------------

* In VirtualBox > Settings > Storage > IDE (Empty)

  - Controller: IDE
  - Choose ISO Image

* In VirtualBox > Start 

 - Select: Install Ubuntu

 - Don't Select: Download updates while installing
 - Gave username and password

 - Select: Continue
 - Select: Install Now
 - Where are you? UK
 - Detect Keyboard
 - Continue
 - Yes to password protect login
 - When installation is complete:

   + Select: Restart Now (it restarts and removes the CD automatically)
   + It asks whether the CD is removed - ensure CD is removed: Press Enter
   + It should restart

**Snapshot taken at this point**


Install Operating System - CentOS
---------------------------------

* In VirtualBox > Settings > Storage > IDE (Empty)

  - Controller: IDE
  - Choose ISO Image

* In VirtualBox > Start 

 - Select: Install CentOS 7
 - Select: English (United Kingdom) - Continue
 - Select the Partitioning
 - I will configure partitioning
 - Done

* New CentOS mount points (there won't be any existing ones). Create:

 - Swap partition (double the RAM we allocated for the VM - i.e. 12GB)
 - Root partition - 20GB (I think this is where software goes)
 - Home partition - the rest of the space i.e. 18GB (for files)
 - For some reason when I entered these numbers, I got 11.18GB, 18.63GB and 16.76GB for swap, root and home respectively. Although the centos volume said 4096kB free. They are all standard partitions.
 - Press Done
 - Accept changes

* Software installation

 - GNOME Applications
 - Internet Applications
 - Legacy X Window System Compatibility
 - Office Suite and Productivity
 - Compatibility Libraries
 - Development Tools
 - Security Tools

* Kdump:

 - Un-enable kdump
 - Done

* Network and hostname:

 - I don't think I have a static IP, so I just turned the ethernet on
 - Done

* Begin installation

* Root Password is set to something

* User is an administrator and has the same password

* Finish configuration

* Reboot the system

* When it rebooted, I selected English and English(UK)

**Snapshot taken at this point**


Change the Settings of the Virtual Machine
------------------------------------------

* In VirtualBox > Settings > Display

  - 128MB Video Memory for increased speed
  - 3D acceleration

* In VirtualBox > Settings > System

  - Processor: change to 4 (needs virtualisation extensions)
  - Execution Cap: 100%

Run Updates for Xubuntu
-----------------------

* In VirtualBox > Start 
* Run updates for Ubuntu

::

    $ sudo apt-get update

Install DKMS
------------

* Install Dynamic Kernel Module Support – otherwise everytime you update the kernel, you'll have to update Guest Additions

::

    $ sudo apt-get install dkms

Install VirtualBox Guest Additions
----------------------------------

* Install Guest Additions, which gives enhanced mouse pointer and video control and installs Vbox Service.

* In the VirtualBox environment: Devices > Install Guest Additions (this opens a window, if not double click the CD mounted on the dekstop)

* Double click autorun.sh and enter user's password (this then installs Guest Additions)

* Right click the CD on the desktop and eject the volume

* Shut down the VM and restart it for Guest Additions to be installed

**Snapshot taken at this point**

Pass Through Folder
===================

* I haven't tried this, but there are loads of tutorials on this, the best one is this:

`Pass Through Folder <https://www.youtube.com/watch?v=1iUafW6g5o8>`_

* Obtain the name of the Ubuntu Box:

::

   sudo nano /etc/hosts

(Could be VMUbuntu)

* In Xubuntu Guest Machine

::

    $ sudo apt-get update
    $ sudo apt-get samba

* In File Manager: Right click folder > Share Folder (will install Samba) 

* Restart session

* In Ubuntu Guest Machine

::

    $ sudo gedit /etc/samba/smb.conf &

* In that file add a line:

::

    [global]

    workgroup = WORKGROUP
    force user = username (this is the username in Ubuntu VM)

* Restart Samba:

::

    $ sudo restart smbd

* Right click folder > Share Folder (Allow others to create and delete files).

 - Create share
 - Create a file in Ubuntu to share

* In Windows Search:

::

    \\VMUbuntu\Public

* In Windows Explorer > Network (May need to turn on share and discovery)

How to Fix Host/Guest Copy Paste Issues
=======================================

* In Virtual Box:

  - Settings > General > Advanced

     + Shared Clipboard: Bidirectional
     + Drag n' Drop: Bidirectional

* If the version of Virtual Box is not up-to-date with respect to the Xubuntu version, this can sometimes not work, if not just update Virtual Box.

Snapshot Management
====================

* Useful for trial software installation. Remember: just running an operating system creates changes. We can't restore while the machine is powered on.

* In VirtualBox > Snapshots:

  - Right click > Take snapshop
  - Right click > Restore snapshot 
  - Right click > Delete snapshot (to remove large files)

How to Clone VMs
================

* How to give other people the VM? **Clone it** Can't give just snapshot

* Clone snapshot: Right click > Clone

**If you are going to be using the clone on the same network as the virtual machine check reinitilise the MAC address of all network cards – otherwise you will have multiple machines with the same MAC address and this will confuse the network.**

* Select: Full clone – so other people can use it
* Select: Clone everything

Can now give other people the .vdi file (the harddisk) and the .vbox. They can launch this directly in their copy of virtual box.

Adjust Network Adapter
======================

* Use: NAT (allows Guest to use Internet)

* Adapter Type: PCnet-Fast III (Am79C973)

* Can randomise MAC address. 08002731A607 (if other VMs are somehow using the same address)

Types of Storage
================

* In VirtualBox: Settings > Storage
* Ubuntu uses SATA, CD uses IDE.

How to Copy vdi files
=====================

* File > Virtual Media Manager
* Differencing files are roughly equivalent to snapshots
* .vdi files have unique identifiers – can't just copy and paste
* **If you want to copy a snapshot – you must use clone (as above)**
* If you want to copy a VM (that isn't a snapshot), you can use: File > Virtual Media Manager > Copy > Next > Next > Next > Give it a name > Copy
* To test the copy, you must add the .vdi file as a **harddrive to a virtual machine** (not just in the virtual machine manager)

Software and Customisation
==========================

* This is a list of software I've installed so far on the VM, plus additional customisation.

::

	# Guest Additions Installed

	ran autorun.sh on Guest Additions CD

	#OpenFOAM installed

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




