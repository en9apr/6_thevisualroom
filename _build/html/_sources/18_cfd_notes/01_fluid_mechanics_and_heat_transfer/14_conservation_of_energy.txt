======================
Conservation of Energy
======================

.. contents::
   :local:

.. highlight:: latex

What is the conservation of total energy in integral form?
==========================================================

* Rate of change of total energy:

.. math::

    {\partial \over {\partial t}} \int_V \rho E dV
    
* Flux of total energy:

.. math::

    \int_S \rho E (\vec{u} \cdot \vec{n}) dS
    
* Work done **by** fluid by surface forces (work is done in the positive direction of :math:`\vec{n}`):

.. math::

    \int_S \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{u} \cdot \vec{n} dS
    
* Heat added **to** fluid by surface forces (heat is added in the negative direction of :math:`\vec{n}`):

.. math::

   -\int_S \vec{q} \cdot \vec{n} dS

* Conservation of energy (neglecting gravity):

.. math::

    {\partial \over {\partial t}} \int_V \rho E dV +
    \int_S \rho E (\vec{u} \cdot \vec{n}) dS = 
    - \int_S \vec{q} \cdot \vec{n} dS +
    \int_S \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{u} \cdot \vec{n} dS
    
What is the conservation of total energy in differential form?
==============================================================

* Gauss divergence theorem:

.. math::

    {\partial \over {\partial t}} \int_V \rho E dV +
    \int_V \nabla \cdot (\rho E \vec{u}) dV = 
    - \int_V \nabla \cdot \vec{q} dV +
    \int_V \nabla \cdot (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) dV

* Volume fixed and arbitary gives:

.. math::

    {\partial \over {\partial t}} (\rho E) +
    \nabla \cdot (\rho E \vec{u}) = 
    - \nabla \cdot \vec{q} +
    \nabla \cdot (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n})
    
* Hence:

.. math::

    {D \over {D t}} (\rho E) = 
    - \nabla \cdot \vec{q} +
    \nabla \cdot (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n})
    
What is the kinetic energy equation in integral form?
=====================================================

* Rate of change of kinetic energy:

.. math::

    {\partial \over {\partial t}} \int_V \rho E_k dV
    
* Flux of kinetic energy:

.. math::

    \int_S \rho E_k (\vec{u} \cdot \vec{n}) dS
    
* Work done **by** fluid by surface forces (work is done in the positive direction of :math:`\vec{n}`):

.. math::

    \int_S \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{u} \cdot \vec{n} dS
    
* Work done **by** fluid by volumetric forces (work is done in the positive direction of :math:`\vec{n}`):

.. math::

   \int_V p \nabla \cdot \vec{u} dV

* Work done **on** fluid by friction (work is done in the negative direction of :math:`\vec{n}`):

.. math::

    - \int_V \Phi dV

* Kinetic energy equation:

.. math::

    {\partial \over {\partial t}} \int_V \rho E_k dV +
    \int_S \rho E_k (\vec{u} \cdot \vec{n}) dS = 
    \int_S \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{u} \cdot \vec{n} dS +
    \int_V p \nabla \cdot \vec{u} dV -
    \int_V \Phi dV
    
What is the kinetic energy equation in differential form?
=========================================================

* Gauss divergence theorem:

.. math::

    {\partial \over {\partial t}} \int_V \rho E_k dV +
    \int_V \nabla \cdot (\rho E_k \vec{u}) dV = 
    \int_V \nabla \cdot (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) dV +
    \int_V p \nabla \cdot \vec{u} dV -
    \int_V \Phi dV
    
* Volume fixed and arbitary gives:

.. math::

    {\partial \over {\partial t}} (\rho E_k) +
    \nabla \cdot (\rho E_k \vec{u}) = 
    \nabla \cdot (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) +
    p \nabla \cdot \vec{u} -
    \Phi 

* Hence:

.. math::

    {D \over {D t}} (\rho E_k) = 
    \nabla \cdot (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) +
    p \nabla \cdot \vec{u} -
    \Phi

How is the internal energy equation derived from the total energy equation and the kinetic energy equation?
===========================================================================================================

* By definition:

.. math::

    {{D} \over {Dt}} (\rho e) = {{D} \over {Dt}} (\rho E) - {{D } \over {Dt}} (\rho E_k)
    
* By substitution:

.. math::

    {{\partial } \over {\partial t}} (\rho e) + \nabla \cdot (\rho e \vec{u}) = 
    - \nabla \cdot \vec{q} +
    \nabla \cdot (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) - 
    (\nabla \cdot (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) +
    p \nabla \cdot \vec{u} -
    \Phi)
    
* By cancellation of shear stress:

.. math::

    {{\partial } \over {\partial t}} (\rho e) + \nabla \cdot (\rho e \vec{u}) = 
    - \nabla \cdot \vec{q} +
    - p \nabla \cdot \vec{u} +
    \Phi
    
* where: :math:`\Phi` = viscous dissipation function

What is the viscous dissipation function?
=========================================

* This is the viscous work put into fluid element deformation and is irreversible.

* It represents the irreversible conversion of mechanical energy to thermal energy through the action of viscosity.

.. math::

    \Phi = 2 \mu \left( e_{ij} - {1 \over 3} \delta_{ij} {{\partial u_j} \over {\partial x_j}} \right)^2

* By expansion:    
    
.. math::

    \Phi = 2 \mu \left[ \left( {{\partial u_1} \over {\partial x_1}} \right)^2 +
                        \left( {{\partial u_2} \over {\partial x_2}} \right)^2 +
                        \left( {{\partial u_3} \over {\partial x_3}} \right)^2 +
                        {1 \over 2} \left( {{\partial u_2} \over {\partial x_1}} + {{\partial u_1} \over {\partial x_2}} \right)^2 +
                        {1 \over 2} \left( {{\partial u_3} \over {\partial x_2}} + {{\partial u_2} \over {\partial x_3}} \right)^2 +
                        {1 \over 2} \left( {{\partial u_1} \over {\partial x_3}} + {{\partial u_3} \over {\partial x_1}} \right)^2
                        \right]
                        
* where:

:math:`\Phi` = viscous dissipation function = :math:`{1 \over 2} \mu \gamma^2`

:math:`e_{ij}` = strain rate tensor = :math:`{1 \over 2} \left( {{\partial u_j} \over {\partial x_i}} + {{\partial u_i} \over {\partial x_j}}   \right)`

:math:`{1 \over 3} \delta_{ij} {{\partial u_j} \over {\partial x_j}}` = compressibility term

What is Fourier's Law?
======================

.. math::

    \vec{q} = -k \nabla T
    
where:

:math:`q` = heat flux per unit area :math:`(W/m^2)`

:math:`k` = thermal conductivity :math:`(W/mK)`

:math:`\nabla T` = temperature gradient :math:`K/m`

* Used in the Navier-Stokes equations to reduce the vector :math:`vec{q}` of three unknowns to one unknown, temperature :math:`T`

* This applies to fluids and solids

What is thermal conductivity?
=============================

* Thermal conductivity is the property of a material to conduct heat

* For water it's a fifth order polynomial

How is the heat transport equation derived from the internal energy equation for an inviscid incompressible flow?
=================================================================================================================

* From the internal energy equation:

.. math::

    {{\partial } \over {\partial t}} (\rho e) + \nabla \cdot (\rho e \vec{u}) = 
    - \nabla \cdot \vec{q} +
    - p \nabla \cdot \vec{u} +
    \Phi
    
where:

.. math::

    e = C_v T
    
and

.. math::

    \vec{q} = -k \nabla T
    
* For an incompressible, inviscid flow: :math:`\nabla \vec{u} = 0`, :math:`\Phi = 0`

.. math::

    {D \over {Dt}} (T) = {k \over {\rho C_v}} \nabla^2 T = \alpha \nabla^2 T
    
* Thermal diffusivity :math:`\alpha` plays the role of viscosity for temperature:

.. math::

    \alpha = {k \over {\rho C_v}}
    
or:

.. math::

    \alpha = {k \over {\rho C_p}}

(as :math:`C_v = C_p` for incompressible flow)

How is the entropy production equation derived from the internal energy equation for incompressible flow?
=========================================================================================================

* From the internal energy equation:

.. math::

    {{\partial } \over {\partial t}} (\rho e) + \nabla \cdot (\rho e \vec{u}) = 
    - \nabla \cdot \vec{q} +
    - p \nabla \cdot \vec{u} +
    \Phi
    
* For a reversible process:

.. math::

    T ds = de + pdv
    
* For incompressible flow:

.. math::

    dv = 0
    
.. math::

    \rho = \text{constant}
    
.. math::

    \nabla \cdot \vec{u} = 0
    
* Hence:

.. math::

    \rho {{De} \over {Dt}} = \rho T {{Ds} \over {Dt}} = -\nabla \cdot \vec{q} + \Phi
    
* Hence:

.. math::

    \rho {{Ds} \over {Dt}} = - {{\nabla \cdot \vec{q}} \over T} + {\Phi \over T}
    
What is the cause of entropy production?
========================================

* Heat conduction (heat added **to** the fluid)

* Viscous dissipation (work done **by** the fluid)