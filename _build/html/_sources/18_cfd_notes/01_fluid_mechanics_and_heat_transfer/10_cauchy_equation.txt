========================
Cauchy Momentum Equation
========================

.. contents::
   :local:

.. highlight:: latex

What is the derivation of the Cauchy Equation?
==============================================

* Temporal change in momentum:

.. math::

    {\partial \over {\partial t}} \int_V \rho \vec{u} dV
    
* Spatial change in momentum:

.. math::

    \int_S \rho \vec{u} (\vec{u} \cdot \vec{n}) dS
    
* Surface forces:

.. math::

    \int_S (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) dS
    
* Body forces:

.. math::

    \int_V \rho \vec{g} dV
    
* Newton's 2nd Law:

.. math::

    \text{Temporal change in momentum + Spatial change in momentum = Surface forces + Body forces}
    
.. math::

    {\partial \over {\partial t}} \int_V \rho \vec{u} dV +
    \int_S \rho \vec{u} (\vec{u} \cdot \vec{n}) dS =
    \int_S (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) dS +
    \int_V \rho \vec{g} dV
    
What is the derivation of the conservation of momentum in differential form?
============================================================================

* Gauss divergence theorem:

.. math::

    \int_S \rho \vec{u} (\vec{u} \cdot \vec{n}) = \int_V \nabla \cdot (\rho \vec{u} \otimes \vec{u}) dV
    
.. math::

    \int_S (\overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} \cdot \vec{n}) dS = \int_V \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} dV


* From integral form:

.. math::

    {\partial \over {\partial t}} \int_V \rho \vec{u} dV +
    \int_V \nabla \cdot (\rho \vec{u} \otimes \vec{u}) dV =
    \int_V \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} dV +
    \int_V \rho \vec{g} dV
    
* Volume is fixed and arbitary:

.. math::

    {\partial \over {\partial t}} (\rho \vec{u}) +
    \nabla \cdot (\rho \vec{u} \otimes \vec{u}) =
    \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} +
    \rho \vec{g}
    
* In Einstein notation:

.. math::

    {\partial \over {\partial t}} (\rho u_j) +
    {\partial \over {\partial x_i}} (\rho u_j u_i) =
    {\partial \over {\partial x_i}} {\sigma_{ij}} +
    \rho g_j
    
What is the non-conservative form of the momentum equation (Cauchy Equation) for incompressible flow?
=====================================================================================================

* Differential form:

.. math::

    {\partial \over {\partial t}} (\rho \vec{u}) +
    \nabla \cdot (\rho \vec{u} \otimes \vec{u}) =
    \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} +
    \rho \vec{g}
    
* Differentiation:

.. math::

    \vec{u} {{\partial \rho} \over {\partial t}} +
    \rho {{\partial \vec{u}} \over {\partial t}} +
    (\vec{u} \cdot \nabla) \rho \vec{u} +
    (\nabla \cdot \vec{u}) \rho \vec{u} =
    \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} +
    \rho \vec{g}
    
* Incompressible :math:`\longrightarrow` :math:`\partial \rho / \partial t = 0`, :math:`\rho = constant` and :math:`\nabla \cdot \vec{u}=0`:

.. math::

    \rho {{\partial \vec{u}} \over {\partial t}}+ (\vec{u} \cdot \nabla) \rho \vec{u} = \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} + \rho \vec{g}

.. math::

    \rho \left( {{\partial \vec{u}} \over {\partial t}}+ (\vec{u} \cdot \nabla) \vec{u} \right) = \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} + \rho \vec{g}  
    
    
What is represented by the Cauchy Equation?
===========================================

* From non-conservative, differential form:

.. math::

    \rho \left( {{\partial \vec{u}} \over {\partial t}}+ (\vec{u} \cdot \nabla) \vec{u} \right) = \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} + \rho \vec{g} 
    
* :math:`\rho {{\partial \vec{u}} \over {\partial t}}` = Temporal change in momentum/unit volume at a fixed point
* :math:`\rho (\vec{u} \cdot \nabla) \vec{u}` = Spatial change in momentum/unit volume at a fixed point
* :math:`\nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma}` = Total stress force/unit volume
* :math:`\rho \vec{g}` = Gravity force/unit volume






