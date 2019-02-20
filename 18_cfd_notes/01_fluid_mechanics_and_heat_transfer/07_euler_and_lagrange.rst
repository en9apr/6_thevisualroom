==================
Euler and Lagrange
==================

.. contents::
   :local:

.. highlight:: latex

What is the difference between the Euler and Lagrange approach?
===============================================================

* Lagrange:
    + Moving frame of reference
    + Elements are matter
    + Results in pathlines
    
* Euler:
    + Stationary frame of reference
    + Elements are spatial
    + Results in streamlines
    
What are the advantages and disadvantages between the Euler and Lagrange approach?
==================================================================================

* Lagrange:
    + Advantage: More suitable than Euler for the dynamics of discrete particles in a fluid e.g. flow visualisation.
    + Disadvantage: Computationally expensive to keep track of large numbers of particles in a flow field

* Euler:
    + Advantage: More suitable than Lagrange for the dynamics of a fluid flow field, e.g. CFD
    + Disadvantage: Time step limited by grid due to stability or accuracy

What is the derivation of the link between the Euler and Lagrange mass conservation?
====================================================================================

* Eulerian mass conservation:

.. math::

    {{\partial \rho} \over {\partial t}} + \nabla \cdot (\rho \vec{u}) = 
    {{\partial \rho} \over {\partial t}} + (\vec{u} \cdot \nabla) \rho  + (\nabla \cdot \vec{u})\rho 
    = 0

* Lagrangian derivative (or Total/Material/Substantive derivative):

.. math::

    {{D \rho} \over {D t}} = {{\partial \rho} \over {\partial t}} + \nabla \cdot (\rho \vec{u}) =
    {{\partial \rho} \over {\partial t}} + (\vec{u} \cdot \nabla) \rho  + (\nabla \cdot \vec{u})\rho 
    
* Lagrangian derivative applies to incompressible flow :math:`\longrightarrow \nabla \cdot \vec{u} = 0`

.. math::

    {{D \rho} \over {D t}} = {{\partial \rho} \over {\partial t}} + (\vec{u} \cdot \nabla) \rho

* Substituting Lagrangian derivative into expanded Eulerian mass conservation:

.. math::

    {{D \rho} \over {D t}} + (\vec{u} \cdot \nabla) \rho = 0
    
