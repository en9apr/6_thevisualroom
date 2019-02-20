============================
Vorticity Transport Equation
============================

.. contents::
   :local:

.. highlight:: latex

What does irrotational mean?
============================

.. math::

    \vec{\Omega} = \vec{0} \longrightarrow \nabla \times \vec{u} = \vec{0}
    
What is an irrotational vortex?
===============================

Particles don't rotate about their own axes, but about the axis of the vortex.

Does irrotational and inviscid, irrotational and viscous, vortical and inviscid, vortical and viscous exist?
============================================================================================================

.. list-table::
   :header-rows: 1
   :widths: 50 30

   * - Type of flow
     - Existing
   * - Irrotational and inviscid
     - Yes
   * - Irrotational and viscous
     - No
   * - Vortical and inviscid
     - Yes
   * - Vortical and viscous
     - Yes

What is circulation?
====================

The line integral around the contour of the tangential component of velocity :math:`\vec{u}`. Via Stokes theorem:

.. math:: \Gamma =  \oint_C \vec{u} \cdot d \vec{C} = \int_S \nabla \times \vec{u} \cdot d \vec{S} = \int_S \vec{\Omega} \cdot d \vec{S}

What is the curl of the gradient of a scalar?
=============================================

The curl of the gradient of a scalar is **zero**:

.. math::

    \nabla \times (\nabla \phi) = \vec{0}
  
Or:

.. math::

    \operatorname{curl} \left( \operatorname{grad} \phi \right) = \vec{0}
  
where :math:`\phi` is a scalar

* If the gradient points in the direction of increasing :math:`\phi`, then
* Travelling in the direction of increasing :math:`\phi` does not rotate.

What is the vector identity for :math:`\vec{u} \cdot \nabla \vec{u}`?
=====================================================================

.. math::

    \vec{u} \cdot \nabla \vec{u} = \nabla \left( {1 \over 2} \vec{u} \cdot \vec{u} \right) - \vec{u} \times \left( \nabla \times \vec{u} \right)
    
What is the curl of :math:`\vec{u} \times \vec{\Omega}` for an incompressible fluid?
====================================================================================

.. math::

    \nabla \times (\vec{u} \times \vec{\Omega}) = \vec{\Omega} \cdot \nabla \vec{u} - \vec{u} \cdot \nabla \vec{\Omega}
    
where: :math:`\nabla \cdot \vec{\Omega} = \nabla \cdot \vec{u} = 0`

How can the non-conservative form of the Navier-Stokes momentum equation use conservative body forces?
======================================================================================================

Navier-Stokes:

.. math::

    {{\partial \vec{u}} \over {\partial t}} + \vec{u} \cdot \nabla \vec{u} = - {1 \over \rho} \nabla p + \nu \nabla^2 \vec{u} + \vec{g}
    

Conservative body forces :math:`\longrightarrow \vec{g} = -\nabla \Phi` 

where: :math:`\Phi = ` gravitational potential function

.. math::

    {{\partial \vec{u}} \over {\partial t}} + \vec{u} \cdot \nabla \vec{u} = - {1 \over \rho} \nabla p + \nu \nabla^2 \vec{u} -\nabla \Phi
    
How can the Navier-Stokes momentum equation be arranged to include the potential function?
==========================================================================================

Navier-Stokes (with conservative body forces):

.. math::

    {{\partial \vec{u}} \over {\partial t}} + \vec{u} \cdot \nabla \vec{u} = - {1 \over \rho} \nabla p + \nu \nabla^2 \vec{u} -\nabla \Phi
 
By expansion: 
 
.. math::

    \vec{u} \cdot \nabla \vec{u} = \nabla \left( {1 \over 2} \vec{u} \cdot \vec{u} \right) - \vec{u} \times \nabla \times \vec{u} = \nabla \left( {1 \over 2} \vec{u} \cdot \vec{u} \right) - \vec{u} \times \vec{\Omega}
    
By substitution:

.. math::

    {{\partial \vec{u}} \over {\partial t}} + \nabla( {1 \over 2} \vec{u} \cdot \vec{u}) - \vec{u} \times \vec{\Omega} = - {1 \over \rho} \nabla p + \nu \nabla^2 \vec{u} -\nabla \Phi

Re-arranging:    
    
.. math::

    {{\partial \vec{u}} \over {\partial t}} - \vec{u} \times \vec{\Omega} = - \nabla \left( {1 \over \rho} p  + {1 \over 2} \vec{u}^2 + \Phi \right) + \nu \nabla^2 \vec{u}

By re-definition:    
    
.. math::

    {{\partial \vec{u}} \over {\partial t}} - \vec{u} \times \vec{\Omega} = - \nabla \Pi + \nu \nabla^2 \vec{u}
     
How can the Navier-Stokes equations using the potential function be expressed as the vorticity transport equation?
==================================================================================================================

Navier-Stokes (with potential function):

.. math::

    {{\partial \vec{u}} \over {\partial t}} - \vec{u} \times \vec{\Omega} = - \nabla \Pi + \nu \nabla^2 \vec{u}
    
Taking the curl of both sides:

.. math::

    \nabla \times \left( {{\partial \vec{u}} \over {\partial t}} - \vec{u} \times \vec{\Omega} \right) = \nabla \times \left( - \nabla \Pi + \nu \nabla^2 \vec{u} \right)

First part:

.. math::

    \nabla \times {{\partial \vec{u}} \over {\partial t}} = {{\partial \vec{\Omega}} \over {\partial t}}
    
Second part:

.. math::

    \nabla \times (- \vec{u} \times \vec{\Omega}) = 
    -\nabla \times \vec{u} \times \vec{\Omega} = 
    - \left( (\vec{\Omega} \cdot \nabla) \vec{u} - (\vec{u} \cdot \nabla) \vec{\Omega}  \right) =
    -(\Omega \cdot \nabla) \vec{u} + (\vec{u} \cdot \nabla) \vec{\Omega}
    
Third part:

.. math::

    \nabla \times (\nabla \Pi) = 0
    
Fourth part:

.. math::

    \nabla \times ( \nu \nabla^2 \vec{u} ) = \nu \nabla^2 \vec{\Omega}
    
Hence:

.. math::

    {{\partial \vec{\Omega}} \over {\partial t}} + (\vec{u} \cdot \nabla) \vec{\Omega} = (\vec{\Omega} \cdot \nabla) \vec{u} + \nu \nabla^2 \vec{\Omega}
    
What is the physical meaning of the terms in the vorticity transport equation?
==============================================================================

.. math::

    {{\partial \vec{\Omega}} \over {\partial t}} + (\vec{u} \cdot \nabla) \vec{\Omega} = (\vec{\Omega} \cdot \nabla) \vec{u} + \nu \nabla^2 \vec{\Omega}
    
* Temporal change of vorticity / unit volume:

.. math:: {{\partial \vec{\Omega}} \over {\partial t}}

* Spatial change of vorticity / unit volume:

.. math:: (\vec{u} \cdot \nabla) \vec{\Omega}

* Vortex stretching term:

.. math:: (\vec{\Omega} \cdot \nabla) \vec{u}

* Diffusion of vorticity / unit volume:

.. math:: \nu \nabla^2 \vec{\Omega}

What is the vortex stretching phenomenon?
=========================================

* The vortex stretching term is :math:`(\vec{\Omega} \cdot \nabla) \vec{u}`

* Conservation of angular momentum:

.. math::

    {{D(I \vec{\Omega})} \over {Dt}} = I {{D \vec{\Omega}} \over {Dt}} + \vec{\Omega} {{DI} \over {Dt}}
    
* For constant angular momentum:

.. math::

    {{D(I \vec{\Omega})} \over {Dt}} = 0
    
* Hence:

.. math::

    {{D \vec{\Omega}} \over {Dt}} = - {{\vec{\Omega}} \over {I}} {{DI} \over {Dt}}
    
* Vortex stretching is an increase in vorticity caused by reduced momentum intertia of vortex elements in 3D

What is the vorticity transport equation in 2D?
===============================================

* In 2D:

.. math::

    (\vec{\Omega} \cdot \nabla) \vec{u} = \left(  \Omega_x {\partial \over {\partial x}} +
    \Omega_y {\partial \over {\partial y}} +
    \Omega_z {\partial \over {\partial z}}
    \right) \vec{u} = 0
    
as :math:`\Omega_x = \Omega_y = {\partial \over {\partial z}} = 0`

* Hence:

.. math::

    {{\partial \vec{\Omega}} \over {\partial t}} + (\vec{u} \cdot \nabla) \vec{\Omega} = \nu \nabla^2 \vec{\Omega}
    
* For :math:`\Omega_z`:

.. math::

    {{\partial \Omega_z} \over {\partial t}} + (\vec{u} \cdot \nabla) \Omega_z = \nu \nabla^2 \Omega_z
    
* By expansion:

.. math::

    {{\partial \Omega_z} \over {\partial t}} + u {{\partial \Omega_z} \over {\partial x}} + v {{\partial \Omega_z} \over {\partial y}}  = \nu \left( {{\partial^2 \Omega_z} \over {\partial x^2}} + {{\partial^2 \Omega_z} \over {\partial y^2}} \right)
   
* Also:

.. math::

    \Omega_x = \Omega_y = 0
    
How are the terms in the 3D vorticity transport equation expanded?
==================================================================

* :math:`\Omega_x` by expansion:

.. math::

    {{\partial \Omega_x} \over {\partial t}} + 
    u {{\partial \Omega_x} \over {\partial x}} + 
    v {{\partial \Omega_x} \over {\partial y}} +
    w {{\partial \Omega_x} \over {\partial z}}
    = \Omega_x {{\partial u} \over {\partial x}} +
    \Omega_y {{\partial u} \over {\partial y}} +
    \Omega_z {{\partial u} \over {\partial z}} +
    \nu \left( 
    {{\partial^2 \Omega_x} \over {\partial x^2}} + 
    {{\partial^2 \Omega_x} \over {\partial y^2}} +
    {{\partial^2 \Omega_x} \over {\partial z^2}}
    \right)

* :math:`\Omega_y` by expansion:

.. math::

    {{\partial \Omega_y} \over {\partial t}} + 
    u {{\partial \Omega_y} \over {\partial x}} + 
    v {{\partial \Omega_y} \over {\partial y}} +
    w {{\partial \Omega_y} \over {\partial z}}
    = \Omega_x {{\partial v} \over {\partial x}} +
    \Omega_y {{\partial v} \over {\partial y}} +
    \Omega_z {{\partial v} \over {\partial z}} +
    \nu \left( 
    {{\partial^2 \Omega_y} \over {\partial x^2}} + 
    {{\partial^2 \Omega_y} \over {\partial y^2}} +
    {{\partial^2 \Omega_y} \over {\partial z^2}}
    \right)  
    
* :math:`\Omega_z` by expansion:

.. math::

    {{\partial \Omega_z} \over {\partial t}} + 
    u {{\partial \Omega_z} \over {\partial x}} + 
    v {{\partial \Omega_z} \over {\partial y}} +
    w {{\partial \Omega_z} \over {\partial z}}
    = \Omega_x {{\partial w} \over {\partial x}} +
    \Omega_y {{\partial w} \over {\partial y}} +
    \Omega_z {{\partial w} \over {\partial z}} +
    \nu \left( 
    {{\partial^2 \Omega_z} \over {\partial x^2}} + 
    {{\partial^2 \Omega_z} \over {\partial y^2}} +
    {{\partial^2 \Omega_z} \over {\partial z^2}}
    \right)  