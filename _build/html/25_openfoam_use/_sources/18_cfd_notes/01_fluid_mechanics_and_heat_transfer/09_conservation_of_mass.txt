====================
Conservation of Mass
====================

.. contents::
   :local:

.. highlight:: latex
    
What is the difference between compressible and incompressible flows?
=====================================================================

* Incompressible flow:
    + :math:`p=f(u)`
    + :math:`\partial p`, :math:`\partial \rho` and :math:`\partial T` (small changes)
* Compressible flow:
    + :math:`p=f(\rho, T)`
    + :math:`\Delta p`, :math:`\Delta \rho` and :math:`\Delta T` (large changes)           
    
What is the derivation of the conservation of mass in integral form?
====================================================================

* Temporal change of mass:

.. math::

    {\partial \over {\partial t}} \int_V \rho dV
    
* Spatial change of mass:

.. math::

    \int_S \rho(\vec{u} \cdot \vec{n})dS

* Conservation of mass:

.. math::

    \text{Temporal change of mass} + \text{Spatial change of mass} = 0
    
.. math::

    {\partial \over {\partial t}} \int_V \rho dV + \int_S \rho(\vec{u} \cdot \vec{n})dS = 0

* If incompressible :math:`\partial \rho / \partial t = 0` and :math:`\rho=constant`:

.. math::

    \int_S (\vec{u} \cdot \vec{n})dS = 0
        
What is the derivation of the conservation of mass in differential form?
========================================================================

* Gauss divergence theorem:

.. math::

    \int_S \rho (\vec{u} \cdot \vec{n})dS = \int_V \nabla \cdot (\rho \vec{u}) dV

* From Integral Form:

.. math::

    {\partial \over {\partial t}} \int_V \rho dV + \int_V \nabla \cdot (\rho \vec{u})dS = 0

* Volume fixed:

.. math::

    \int_V \left( {{\partial \rho} \over {\partial t}} + \nabla \cdot (\rho \vec{u}) \right) dV = 0
    
* Volume is arbitary (this is for compressible flow):

.. math::

    {{\partial \rho} \over {\partial t}} + \nabla \cdot (\rho \vec{u}) = 0
    
* By expansion:

.. math::

    {{\partial \rho} \over {\partial t}} + (\vec{u} \cdot \nabla) \rho  + (\nabla \cdot \vec{u}) \rho  = 0
    
* If incompressible :math:`{{\partial \rho} \over {\partial t}} = 0` and :math:`\nabla \rho = 0` and :math:`\rho = constant`:

.. math::

    \nabla \cdot \vec{u} = 0