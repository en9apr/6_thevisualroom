============================
How to Use a Virtual Machine
============================

.. contents::
   :local:


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
