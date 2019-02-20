==========
Kinematics
==========

.. contents::
   :local:

.. highlight:: latex

What is the main difference between kinematics and dynamics?
============================================================

* Kinematics:

    - Quantities involving space and time only, e.g.
        + Position
        + Velocity
        + Acceleration
        + Deformation
        + Rotation
    - i.e. the geometry of fluid motion

* Dynamics:

    - Refers to the stresses and forces that cause fluid motion
    - i.e. the conservation laws
    
What are the different types of motion exhibited by fluid elements?
===================================================================

* Translation :math:`\vec{u}`
* Rotation :math:`{1 \over 2} (\nabla \times \vec{u}) \times d\vec{r}`
* Angular deformation or linear deformation :math:`\vec{s} \cdot d\vec{r}`

.. math::

    \vec{u} + {1 \over 2} (\nabla \times \vec{u}) \times d\vec{r} + \vec{s} \cdot d \vec{r}
    

How can fluid flow be decomposed into translation, deformation and rotation?
============================================================================

.. math::

    \vec{u} = \vec{u} + 
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{e} \cdot d\vec{c} +
    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\omega} \times d\vec{c}
    = translation + deformation + rotation
    
where:

* :math:`\vec{u}` = velocity vector
* :math:`\overset{\underset{\mathrm{\rightrightarrows}}{}}{e}` = strain rate tensor
* :math:`\vec{c}` = line element vector
* :math:`\overset{\underset{\mathrm{\rightrightarrows}}{}}{\omega}` = vorticity tensor
 
How to decompose the velocity gradient tensor into symmetric and anti-symmetric parts?
======================================================================================

.. math::

    {{\partial u_j} \over {\partial x_i}} = 
    {1 \over 2} \left( {{\partial u_j} \over {\partial x_i}} + {{\partial u_i} \over {\partial x_j}} \right) +
    {1 \over 2} \left( {{\partial u_j} \over {\partial x_i}} - {{\partial u_i} \over {\partial x_j}} \right) = 
    e_{ij} + \omega_{ij} = symmetric \ tensor + antisymmetric \ tensor
    
:math:`{{\partial u_j} \over {\partial x_i}}` appears in the Navier-Stokes equations (shear strain rate tensor).

.. math::

    e_{ij} = {1 \over 2} \left( {{\partial u_j} \over {\partial x_i}} + {{\partial u_i} \over {\partial x_j}} \right)
    
(eventually this becomes the shear stress)

What is the rate of strain tensor or deformation tensor?
========================================================

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{e} =
    {1 \over 2} \left( \nabla \otimes \vec{u} + (\nabla \otimes \vec{u})^T \right) = 
    {1 \over 2} \left( \nabla_i u_j + \nabla_j u_i \right) =
    \begin{bmatrix}
    {\partial u \over \partial x} & {1 \over 2}\left( {\partial v \over \partial x} + {\partial u \over \partial y} \right)  & {1 \over 2}\left( {\partial w \over \partial x} + {\partial u \over \partial z} \right)  \\
    {1 \over 2}\left( {\partial u \over \partial y} + {\partial v \over \partial x} \right) & {\partial v \over \partial y} & {1 \over 2}\left( {\partial w \over \partial y} + {\partial v \over \partial z} \right) \\
    {1 \over 2}\left( {\partial u \over \partial z} + {\partial w \over \partial x} \right) & {1 \over 2}\left( {\partial v \over \partial z} + {\partial w \over \partial y} \right) & {\partial w \over \partial z} 
    \end{bmatrix}

What is the vorticity tensor?
=============================

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{\omega} =
    {1 \over 2} \left( \nabla \otimes \vec{u} - (\nabla \otimes \vec{u})^T \right) = 
    {1 \over 2} \left( \nabla_i u_j - \nabla_j u_i \right) =
    \begin{bmatrix}
    0 & {1 \over 2}\left( {\partial v \over \partial x} - {\partial u \over \partial y} \right)  & {1 \over 2}\left( {\partial w \over \partial x} - {\partial u \over \partial z} \right)  \\
    {1 \over 2}\left( {\partial u \over \partial y} - {\partial v \over \partial x} \right) & 0 & {1 \over 2}\left( {\partial w \over \partial y} - {\partial v \over \partial z} \right) \\
    {1 \over 2}\left( {\partial u \over \partial z} - {\partial w \over \partial x} \right) & {1 \over 2}\left( {\partial v \over \partial z} - {\partial w \over \partial y} \right) & 0 
    \end{bmatrix}

What is represented by the symmetric and anti-symmetric parts of the velocity gradient tensor?
==============================================================================================

* symmetric part :math:`\longrightarrow` angular deformation :math:`\longrightarrow` shear
* anti-symmetric part :math:`\longrightarrow` rotation :math:`\longrightarrow` vorticity
