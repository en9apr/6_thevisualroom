============================
Integral and Derivative Form
============================

.. contents::
   :local:

.. highlight:: latex
    
What are the conditions of applicability of the integral and differential forms of the conservation laws?
=========================================================================================================

* Integral form :math:`\longrightarrow` same as limits of continuum model:

    + :math:`Kn \lt \lt 1`
    + No mass-energy conversion

* Differential form :math:`\longrightarrow` same as limits of Gauss-Divergence theorem:

    + :math:`\phi` must be continuous (no shock waves or free surfaces etc)
    + :math:`\partial \phi / \partial x` must exist
    + :math:`\partial \phi / \partial x` must be continuous
 
What is the Reynolds Transport Theorem?
=======================================

.. math::

    {D \over {Dt}} \int_V \phi dV = 
    {\partial \over {\partial t}} \int_V \phi dV + \int_S \phi (\vec{u} \cdot \vec{n}) dS

* where: :math:`\phi` is a scalar or vector and is a conserved quantity

What is the meaning of the Reynolds Transport Theorem?
======================================================

.. math::

    \mathbf{Temporal} \ change \ of \ \phi \ following \ a \ \mathbf{moving} \ volume = 
    \mathbf{Temporal} \ change \ of \ \phi \ in \ a \ \mathbf{fixed} \ volume +
    \mathbf{Spatial} \ change \ of \ \phi \ through \ \mathbf{fixed} \ surfaces
    
What are some examples of the Reynolds Transport Theorem?
=========================================================

* Mass conservation :math:`\phi = \rho` (mass per unit volume):

.. math::

    {D \over {Dt}} \int_V \rho dV = 
    {\partial \over {\partial t}} \int_V \rho dV + \int_S \rho (\vec{u} \cdot \vec{n}) dS
    
* Momentum conservation :math:`\phi = \rho \vec{u}` (momentum per unit volume):

.. math::

    {D \over {Dt}} \int_V \rho \vec{u} dV = 
    {\partial \over {\partial t}} \int_V \rho \vec{u} dV + \int_S \rho \vec{u} (\vec{u} \cdot \vec{n}) dS
     
* Energy conservation :math:`\phi = \rho E` (total specific energy per unit volume):

.. math::

    {D \over {Dt}} \int_V \rho E dV = 
    {\partial \over {\partial t}} \int_V \rho E dV + \int_S \rho E (\vec{u} \cdot \vec{n}) dS
        
What is the total derivative?
=============================

* From Reynolds Transport Theorem:

.. math::

    {D \over {Dt}} \int_V \phi dV = 
    {\partial \over {\partial t}} \int_V \phi dV + \int_S \phi (\vec{u} \cdot \vec{n}) dS
    
* Gauss Divergence:

.. math::

    {D \over {Dt}} \int_V \phi dV = 
    {\partial \over {\partial t}} \int_V \phi dV + \int_V \nabla \cdot (\phi \vec{n}) dV

* Equal Volumes:

.. math::

    {{D \phi} \over {Dt}} = {{\partial \phi} \over {\partial t}} + \nabla \cdot (\phi \vec{u})
    
* By expansion:

.. math::

    {{D \phi} \over {Dt}} = {{\partial \phi} \over {\partial t}} + 
    (\vec{u} \cdot \nabla) \phi + (\nabla \cdot \vec{u}) \phi
   
* Incompressible :math:`\longrightarrow \nabla \cdot \vec{u} = 0`:

.. math::

    {{D \phi} \over {Dt}} = {{\partial \phi} \over {\partial t}} + 
    (\vec{u} \cdot \nabla) \phi  = 
    {{\partial \phi} \over {\partial t}} + 
    u_i {{\partial \phi} \over {\partial x_i}}
   
What is the meaning of the total derivative?
============================================

.. math::

    \mathbf{Temporal} \ change \ of \ \phi \ following \ a \ \mathbf{moving} \ point = 
    \mathbf{Temporal} \ change \ of \ \phi \ at \ a \ \mathbf{fixed} \ point +
    \mathbf{Spatial} \ change \ of \ \phi \ at \ a \ \mathbf{fixed} \ point
    
Or

.. math::

    Total \ Acceleration = Unsteady \ Acceleration + Convective/Advective \ Acceleration 
    
    
    
    
    
    
