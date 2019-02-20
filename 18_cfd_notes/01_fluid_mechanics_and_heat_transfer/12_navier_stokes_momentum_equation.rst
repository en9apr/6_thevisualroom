===============================
Navier-Stokes Momentum Equation
===============================

.. contents::
   :local:

.. highlight:: latex

What is the relationship between dynamic and kinematic viscosity and what are their units?
==========================================================================================

.. math::

    \nu = {\mu \over \rho}
    
**Mnemonic:**

Kinematic viscosity :math:`\ll` Dynamic viscosity, as kinematics deals with smaller things than dynamics.

**Units:**

:math:`\nu \ (m^2/s)`

:math:`\mu \ (kg/m/s)`

:math:`\rho \ (kg/m^3)`

What is the effect of the incompressibility constraint on the Navier-Stokes equations?
======================================================================================

1) :math:`\rho = \text{constant} \longrightarrow \nabla \cdot \vec{u} = 0` i.e. velocity field must be divergence free
2) :math:`\rho = \text{constant} \longrightarrow` pressure and density are decoupled, i.e. pressure is no longer a state variable, it is a constant :math:`\longrightarrow` require pressure-velocity coupling

What are the assumptions in the Navier-Stokes equations?
========================================================

1) Continuum hypothesis :math:`\longrightarrow Kn \ll 1` and no mass-energy conversion
2) Form of diffusive fluxes e.g.:
    * Newtownian :math:`\tau_{ij} = \mu \gamma_{ij}`
    * Fourier's Law :math:`\vec{q} = -k \nabla T`
    * Viscosity model :math:`\mu = \mu(T)` e.g. Sutherland's Law
3) Equation of state (but this concerns the solution)
    * Stokes assumption :math:`{2 \over 3} \mu + \lambda = 0`
    * Thermally perfect gas :math:`p=\rho RT`
    * Calorcally perfect gas :math:`e=C_v T`
    
What is the relationship between the total stress tensor, normal stress and deviatoric stress?
==============================================================================================

.. math::

    \text{Total stress tensor = Normal stress tensor + Deviatoric stress tensor}

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} = 
    -p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} + 
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\tau}

(Deviatoric stress tensor is also the shear stress tensor)    
    
where the identity matrix is:

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} = 
    \begin{bmatrix}
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 0 & 1 
    \end{bmatrix}
    
What is the definition of the total stress tensor for viscous and inviscid flow?
================================================================================

Total stress tensor:

* Pressure tensor (normal stress tensor)
* Deviatoric stress tensor (shear stress tensor)

Inviscid:

* Fluid only has normal force:

.. math::
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\tau} = 0

* Hence:

.. math::
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} = -p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I}

Viscous:

* Fluid has shear stress: 

.. math::
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\tau} \ne 0


* Hence:

.. math::
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} = -p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} + \overset{\underset{\mathrm{\rightrightarrows}}{}}{\tau}

How is the total stress tensor defined for the Navier-Stokes momentum equation?
===============================================================================

* Total stress tensor:

.. math::
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma}=-p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} + \overset{\underset{\mathrm{\rightrightarrows}}{}}{\tau}

* Shear stress tensor:

.. math::
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\tau}= \mu \overset{\underset{\mathrm{\rightrightarrows}}{}}{\gamma}
    
* Shear strain rate tensor:

.. math::
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\gamma}= 2 \left( \overset{\underset{\mathrm{\rightrightarrows}}{}}{e} - {1 \over 3} (\nabla \cdot \vec{u}) \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} \right)

* Shear strain rate tensor:

.. math::
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{e}= {1 \over 2} \left( \nabla \otimes \vec{u} + (\nabla \otimes \vec{u})^T \right)

* Hence:

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma}=-p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} + 
    \mu \left( \nabla \otimes \vec{u} + (\nabla \otimes \vec{u})^T - {2 \over 3} (\nabla \cdot \vec{u}) \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} \right)
    
What is the shear strain rate tensor?
=====================================

* In Einstein notation:

.. math::

    e_{ij} = {1 \over 2} \left( {{\partial u_j} \over {\partial x_i}} + {{\partial u_i} \over {\partial x_j}} \right)

.. math::

    e_{ij} = 
    \begin{bmatrix}
    {{\partial u_1} \over {\partial x_1}} &  {1 \over 2} \left( {{\partial u_2} \over {\partial x_1}} + {{\partial u_1} \over {\partial x_2}} \right) & {1 \over 2} \left( {{\partial u_3} \over {\partial x_1}} + {{\partial u_1} \over {\partial x_3}} \right) \\
    {1 \over 2} \left( {{\partial u_1} \over {\partial x_2}} + {{\partial u_2} \over {\partial x_1}} \right) &  {{\partial u_2} \over {\partial x_2}} & {1 \over 2} \left( {{\partial u_3} \over {\partial x_2}} + {{\partial u_2} \over {\partial x_3}} \right) \\
    {1 \over 2} \left( {{\partial u_1} \over {\partial x_3}} + {{\partial u_3} \over {\partial x_1}} \right) &  {1 \over 2} \left( {{\partial u_2} \over {\partial x_3}} + {{\partial u_3} \over {\partial x_2}} \right) & {{\partial u_3} \over {\partial x_3}}
    \end{bmatrix}

* In vector notation:

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{e}= {1 \over 2} \left( \nabla \otimes \vec{u} + (\nabla \otimes \vec{u})^T \right)
    
How is the compressible Navier-Stokes momentum equation derived from the Cauchy equation using the shear strain rate tensor?
============================================================================================================================

* Cauchy:

.. math::

    {D \over {Dt}} (\rho \vec{u}) = \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} + \rho \vec{g}
    
* Total stress tensor:

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} = -p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} +
    \mu \left( \nabla \otimes \vec{u} + (\nabla \otimes \vec{u})^T - {2 \over 3} (\nabla \cdot \vec{u}) \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} \right)
    
* Taking the divergence of each part:

.. math::

    \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} = \nabla \cdot (-p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I}) +
    \mu \nabla \cdot \left( \nabla \otimes \vec{u} + (\nabla \otimes \vec{u})^T \right) - 
    {2 \over 3} \mu \nabla \cdot ( (\nabla \cdot \vec{u}) \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} )

* First part:

.. math::

    \nabla \cdot (-p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I}) = -(\nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{I})p - 
    (\overset{\underset{\mathrm{\rightrightarrows}}{}}{I} \cdot \nabla)p = -\nabla p
    

* Second part:

.. math::

    \mu \nabla \cdot \left( \nabla \otimes \vec{u} + (\nabla \otimes \vec{u})^T \right) = 
    \mu \left(\nabla^2 \vec{u} + \nabla (\nabla \cdot \vec{u}) \right) = 
    \mu \nabla^2 \vec{u} + \mu \nabla (\nabla \cdot \vec{u})
    
* Third part:

.. math::

    -{2 \over 3} \mu \nabla \cdot ( (\nabla \cdot \vec{u}) \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} ) =
    -{2 \over 3} \mu \left( (\nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{I}) \nabla \cdot \vec{u} +
    (\overset{\underset{\mathrm{\rightrightarrows}}{}}{I} \cdot \nabla) \nabla \cdot \vec{u}   \right) = 
    -{2 \over 3} \mu \nabla (\nabla \cdot \vec{u})

where: :math:`\nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} = 0` and :math:`\overset{\underset{\mathrm{\rightrightarrows}}{}}{I} \cdot \nabla=\nabla`

* Adding first, second and third parts (plus :math:`\rho \vec{g}`):

.. math::

    {D \over {Dt}} (\rho \vec{u}) = -\nabla p + \mu \nabla^2 \vec{u} + \mu \nabla (\nabla \cdot \vec{u}) -{2 \over 3} \mu \nabla (\nabla \cdot \vec{u}) + \rho \vec{g}
    
* Simplfiying:

.. math::

    {D \over {Dt}} (\rho \vec{u}) = -\nabla p + \mu \left( \nabla^2 \vec{u} + {1 \over 3} \nabla (\nabla \cdot \vec{u}) \right) + \rho \vec{g}
    
What are the meaning of the terms in the fully compressible Navier-Stokes equations?
====================================================================================

.. math::

    {\partial \over {\partial t}} (\rho \vec{u}) + \nabla \cdot (\rho \vec{u} \otimes \vec{u} ) = -\nabla p + \mu \left( \nabla^2 \vec{u} + {1 \over 3} \nabla (\nabla \cdot \vec{u}) \right) + \rho \vec{g}

* Unsteady momentum:

.. math::

    {\partial \over {\partial t}} (\rho \vec{u})

* Advective momentum:

.. math::

    \nabla \cdot (\rho \vec{u} \otimes \vec{u} )

* Pressure gradient:

.. math::

    - \nabla p

* Viscous or diffusive term:

.. math::

    \mu \nabla^2 \vec{u}

* Compressibility:

.. math::

    {1 \over 3} \mu \nabla (\nabla \cdot \vec{u})

* Body forces:

.. math::

    \rho \vec{g}
    
Expand the terms in the Navier-Stokes equations assuming :math:`\nabla \cdot \vec{u} = 0`
=========================================================================================

* Navier-Stokes:

.. math::

    {D \over {Dt}} (\rho \vec{u}) = {\partial \over {\partial t}} (\rho \vec{u}) + \nabla \cdot (\rho \vec{u} \otimes \vec{u} ) = -\nabla p + \mu \left( \nabla^2 \vec{u} + {1 \over 3} \nabla (\nabla \cdot \vec{u}) \right) + \rho \vec{g}
    
* :math:`x`-direction as principal:

.. math::

    {\partial \over {\partial t}} (\rho u) +
    {\partial \over {\partial x}} (\rho u^2) +
    {\partial \over {\partial x}} (\rho uv) +
    {\partial \over {\partial x}} (\rho uw) =
    -{{\partial p} \over {\partial x}} +
    \mu \left( {{\partial^2 u} \over {\partial x^2}} + 
    {{\partial^2 u} \over {\partial y^2}} + 
    {{\partial^2 u} \over {\partial z^2}}  \right) +
    \rho g_x
    
* :math:`y`-direction as principal:

.. math::

    {\partial \over {\partial t}} (\rho v) +
    {\partial \over {\partial x}} (\rho vu) +
    {\partial \over {\partial x}} (\rho v^2) +
    {\partial \over {\partial x}} (\rho vw) =
    -{{\partial p} \over {\partial y}} +
    \mu \left( {{\partial^2 v} \over {\partial x^2}} + 
    {{\partial^2 v} \over {\partial y^2}} + 
    {{\partial^2 v} \over {\partial z^2}}  \right) +
    \rho g_y

* :math:`z`-direction as principal:

.. math::

    {\partial \over {\partial t}} (\rho w) +
    {\partial \over {\partial x}} (\rho wu) +
    {\partial \over {\partial x}} (\rho wv) +
    {\partial \over {\partial x}} (\rho x^2) =
    -{{\partial p} \over {\partial z}} +
    \mu \left( {{\partial^2 w} \over {\partial x^2}} + 
    {{\partial^2 w} \over {\partial y^2}} + 
    {{\partial^2 w} \over {\partial z^2}}  \right) +
    \rho g_z
    
How is the compressibility term expanded in the Navier-Stokes equations?
========================================================================

* Compressibility term:

.. math::

    {1 \over 3} \nabla (\nabla \cdot \vec{u})
    
* :math:`x`-direction as principal:

.. math::

    {1 \over 3} {\partial \over {\partial x}} \left(
    {{\partial u} \over {\partial x}} + 
    {{\partial v} \over {\partial y}} + 
    {{\partial w} \over {\partial z}} \right)
    
* :math:`y`-direction as principal:

.. math::

    {1 \over 3} {\partial \over {\partial y}} \left(
    {{\partial u} \over {\partial x}} + 
    {{\partial v} \over {\partial y}} + 
    {{\partial w} \over {\partial z}} \right)    
    
* :math:`z`-direction as principal:

.. math::

    {1 \over 3} {\partial \over {\partial z}} \left(
    {{\partial u} \over {\partial x}} + 
    {{\partial v} \over {\partial y}} + 
    {{\partial w} \over {\partial z}} \right)    
    
What is the skew symmetric form of the total derivative?
========================================================

.. math::

    {D \over {Dt}} (\rho \vec{u}) = {{\partial \vec{u}} \over {\partial t}} + {1 \over 2} \left( (\vec{u} \cdot \nabla) \vec{u} + \nabla \cdot (\vec{u} \otimes \vec{u})  \right)