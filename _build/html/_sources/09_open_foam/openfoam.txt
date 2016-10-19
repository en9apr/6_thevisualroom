=================
OpenFOAM Commands
=================

Here is a list of some commands in OpenFOAM

.. contents::
   :local:

How to build a mesh?
====================

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Command
     - Meaning
   * - ::

           $ blockMesh

     - Creates parametric meshes with grading and curved edges

How to run a case?
==================

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Command
     - Meaning
   * - ::

           $ icoFoam

     - Laminar, incompressible, isothermal solver
   * - ::

           $ pisoFoam

     - Turbulent, incompressible, isothermal solver


How to post-process data?
=========================

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Command
     - Meaning
   * - ::

           $ paraFoam

     - Opens Paraview

