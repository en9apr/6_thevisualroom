=======================
Euler Momentum Equation
=======================

.. contents::
   :local:

.. highlight:: latex

What is the derivation of the Euler Momentum Equation from the Cauchy Equation?
===============================================================================

* Cauchy Equation:

.. math::

    \rho \left( {{\partial \vec{u}} \over {\partial t}}+ (\vec{u} \cdot \nabla) \vec{u} \right) = \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} + \rho \vec{g} 
    
* Inviscid flow :math:`\longrightarrow \overset{\underset{\mathrm{\rightrightarrows}}{}}{\tau}=0`:

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} = 
    \begin{bmatrix}
    \sigma_{11} & 0 & 0 \\
    0 & \sigma_{22} & 0 \\
    0 & 0 & \sigma_{33} 
    \end{bmatrix}
    
* Normal stresses = applied pressure:

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} = -p \delta_{ij} 
    = -p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I} = 
    \begin{bmatrix}
    -p & 0 & 0 \\
    0 & -p & 0 \\
    0 & 0 & -p 
    \end{bmatrix}
    
* Expansion of divergence of shear stress:

.. math::

    \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{\sigma} = 
    \nabla \cdot (-p \overset{\underset{\mathrm{\rightrightarrows}}{}}{I}) =
    -(\nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{I})p 
    -(\overset{\underset{\mathrm{\rightrightarrows}}{}}{I} \cdot \nabla)p = -\nabla p

* Hence:

.. math::

    \rho \left( {{\partial \vec{u}} \over {\partial t}}+ (\vec{u} \cdot \nabla) \vec{u} \right) = 
    -\nabla p + \rho \vec{g} 
    
What is the meaning of the terms in the Euler momentum equation?
================================================================

.. math::

    \rho \left( {{\partial \vec{u}} \over {\partial t}}+ (\vec{u} \cdot \nabla) \vec{u} \right) = 
    -\nabla p + \rho \vec{g} 
    
* :math:`\rho {{\partial \vec{u}} \over {\partial t}}` = Temporal change of momentum/unit volume at a fixed point
* :math:`\rho(\vec{u} \cdot \nabla) \vec{u}` = Spatial change of momentum/unit volume at a fixed point
* :math:`-\nabla p` = Pressure force/unit volume at a fixed point
* :math:`\rho \vec{g}` = Gravity force/unit volume at a fixed point

How are all the terms in the Euler Equations expanded in 3D?
============================================================

* x-direction:

.. math::

    \rho \left( {{\partial {u}} \over {\partial t}}+ {u} {{\partial u} \over {\partial x}} +
    {v} {{\partial u} \over {\partial y}} + {w} {{\partial u} \over {\partial z}} \right)= 
    -{{\partial {p}} \over {\partial x}} + \rho g_x

* y-direction:

.. math::

    \rho \left( {{\partial {v}} \over {\partial t}}+ {u} {{\partial v} \over {\partial x}} +
    {v} {{\partial v} \over {\partial y}} + {w} {{\partial v} \over {\partial z}} \right)= 
    -{{\partial {p}} \over {\partial y}} + \rho g_y

* z-direction:
    
.. math::

    \rho \left( {{\partial {w}} \over {\partial t}}+ {u} {{\partial w} \over {\partial x}} +
    {v} {{\partial w} \over {\partial y}} + {w} {{\partial w} \over {\partial z}} \right)= 
    -{{\partial {p}} \over {\partial z}} + \rho g_z  
    
    
* Can also add source terms to the RHS :math:`S_x`, :math:`S_y` and :math:`S_z` 