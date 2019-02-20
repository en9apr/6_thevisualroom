==========================
Vector and Tensor Calculus
==========================

.. contents::
   :local:

.. highlight:: latex

What can integral and differential equations describe in fluid mechanics?
=========================================================================

Differential equations describe:

* A **point** or particle

Integral equations describe:

* A **control volume** and fluxes across control surfaces


What kind of information can integral and differential equations provide?
=========================================================================

Differential equations are needed for:

* Field information, i.e. **continuous problems**

Integral equations are needed for:

* Average quantities, i.e. **discontinuous problems**

What is the degree of an equation?
==================================

The highest power, i.e. Equation :eq:`power` is 2nd degree, also called **non-linear**.  

.. math:: {{\partial \over {\partial t}} (\rho u)} + {{\partial \over {\partial x}} (\rho u^2)} = -{{\partial p} \over {\partial x}} + \mu {{\partial^2 u} \over {\partial x^2}}
    :label: power

If the highest power is 1, i.e. Equation :eq:`linear` is 1st degree, also called **linear**.  

.. math:: {{\partial u} \over {\partial t}}  + c {{\partial u} \over {\partial x}}=0
    :label: linear

What is the order of an equation?
=================================

The highest derivative, i.e. Equation :eq:`order` is 2nd order.

.. math:: {{\partial u} \over {\partial t}}  = {{\partial^2 u} \over {\partial x^2}}
    :label: order
    
What is the difference between ODEs and PDEs?
=============================================

ODEs:

* An equation relating a dependent variable (say x) to **one independent variable** (say u).
* Uses ordinary derivative notation e.g. :math:`dx / dt = a` and :math:`du / dt  =0`.
* Suitable for describing fluids in a **moving** frame of reference.

PDEs:

* An equation relating a dependent variable (say u) to **more than one independent variable** (say x and t).
* Uses partial derivative notation e.g. :math:`\partial u / \partial t + a (\partial u / \partial x) = 0`.
* Suitable for describing fluids in a **fixed** frame of reference.

What is the difference between index and invariant notation?
============================================================

Index notation:

* **Coordinate system**, e.g. Cartesian, Cylindrical, Spherical (:math:`x`, :math:`y`, :math:`z`) or (:math:`x_1`, :math:`x_2`, :math:`x_3`)

Invariant notation:

* **Coordinate free system**, e.g. vector notation (div, grad, curl)

What is Stokes' Theorem?
========================

Converts line integrals to surface integrals (:math:`\vec{u}` is continuous).

.. math:: \Gamma =  \oint_C \vec{u} \cdot d \vec{C} = \int_S \nabla \times \vec{u} \cdot d \vec{S} = \int_S \vec{\omega} \cdot d \vec{S}

(It's a bit like the Gauss divergence theorem, it goes from a length to an area - Gauss divergence theorem goes from an area to a volume).

What is the gradient of a scalar field?
=======================================

.. math:: 

    \text{grad } u = \nabla u =  
    \begin{bmatrix}
    \partial / {\partial x_1}  \\
    \partial / {\partial x_2} \\
    \partial / {\partial x_3}
    \end{bmatrix} ^T u = 
    {{\partial u} \over {\partial x_1}} \vec{i} +
    {{\partial u} \over {\partial x_2}} \vec{j} +
    {{\partial u} \over {\partial x_3}} \vec{k} =
    {{\partial u} \over {\partial x_j}}\vec{e}_i = 
    {{\partial u_i} \over {\partial x_j}} = 
    {{\partial u} \over {\partial x}} \vec{i} +
    {{\partial u} \over {\partial y}} \vec{j} +
    {{\partial u} \over {\partial z}} \vec{k} 

* :math:`i` = principal direction
* :math:`j` = number of dimensions
* Always points in the direction of increasing u
* In the Navier-Stokes equations, :math:`\nabla \vec{u}` is sometimes used implying :math:`\vec{u}` is a row vector for each principal direction. 

What is the gradient of a vector?
=================================

.. math:: 

    \text{grad } \vec{u} = \nabla \vec{u}  =  
    \begin{bmatrix}
    \partial / {\partial x_1}  \\
    \partial / {\partial x_2} \\
    \partial / {\partial x_3}
    \end{bmatrix} ^T 
    \begin{bmatrix}
    u_1  \\
    u_2 \\
    u_3
    \end{bmatrix} = 
    {\sum_{j=1}^3} \nabla_j u_i = 
    \begin{bmatrix}
    {{\partial u_1} \over {\partial x_1}} +  {{\partial u_1} \over {\partial x_2}} + {{\partial u_1} \over {\partial x_3}}\\
    {{\partial u_2} \over {\partial x_1}} +  {{\partial u_2} \over {\partial x_2}} + {{\partial u_2} \over {\partial x_3}} \\
    {{\partial u_3} \over {\partial x_1}} +  {{\partial u_3} \over {\partial x_2}} + {{\partial u_3} \over {\partial x_3}}
    \end{bmatrix}
    
What is the curl of a vector field?
===================================

In 3D:

.. math::

    \nabla \times \vec{u} = 
    \begin{vmatrix}
    \vec{i} & \vec{j} & \vec{k}  \\
    \partial / {\partial x_1} & \partial / {\partial x_2} & \partial / {\partial x_3} \\
    u_1 & u_2 & u_3
    \end{vmatrix} =
    \left( {{\partial u_3} \over {\partial x_2}}  -  {{\partial u_2} \over {\partial x_3}} \right) \vec{i} + 
    \left( {{\partial u_1} \over {\partial x_3}} -  {{\partial u_3} \over {\partial x_1}} \right) \vec{j} + 
    \left( {{\partial u_2} \over {\partial x_1}} -  {{\partial u_1} \over {\partial x_2}} \right) \vec{k} 
    
In 2D:

.. math::

    \nabla \times \vec{u} = 
    \begin{vmatrix}
    \partial / {\partial x_1} & \partial / {\partial x_2} \\
    u_1 & u_2
    \end{vmatrix} =
    {{\partial u_2} \over {\partial x_1}} -  {{\partial u_1} \over {\partial x_2}}

Method in 3D:
    
* Cover row and column for each component :math:`i`, :math:`j` and :math:`k`
* Take determinant of :math:`2 \times 2` matrix

The physical meaning is the amount to which the vector field rotates.

What is the divergence of the Kronecker delta?
==============================================


.. math::

    \nabla_i \cdot \delta_{ij} = 
    \begin{bmatrix}
    \partial / {\partial x_1} \\ \partial / {\partial x_2} \\ \partial / {\partial x_3} 
    \end{bmatrix}^T
    \begin{bmatrix}
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 0 & 1 
    \end{bmatrix} =
    \begin{bmatrix}
    \partial / {\partial x_1} & \partial / {\partial x_2} & \partial / {\partial x_3} 
    \end{bmatrix}
    
i.e.

.. math::

    \nabla_i \cdot \delta_{ij} = \nabla_j
    
or

.. math::

    {\partial \over {\partial x_i}} {\delta_{ij}} = {\partial \over {\partial x_j}}
    
This operation converts a column vector (:math:`\nabla_i`) into a row vector (:math:`\nabla_j`).

What is the divergence of a vector field?
=========================================

.. math:: 

    \text{div } \vec{u} = \nabla \cdot \vec{u}  =  
    \begin{bmatrix}
    \partial / {\partial x_1}  \\
    \partial / {\partial x_2} \\
    \partial / {\partial x_3}
    \end{bmatrix} ^T 
    \cdot
    \begin{bmatrix}
    u_1  \\
    u_2 \\
    u_3
    \end{bmatrix} = 
    {{\partial u_1} \over {\partial x_1}} +  {{\partial u_2} \over {\partial x_2}} + {{\partial u_3} \over {\partial x_3}} = 
    {{\partial u_i} \over {\partial x_i}}
    
* The extent to which the vector field is a source (+ve) or sink (-ve)
* If :math:`div \vec{u} = 0` the vector field is divergence free

What is the Gauss divergence theorem?
=====================================

.. math::

    \int_S \phi (\vec{u} \cdot \vec{n}) dS = \int_V \nabla \cdot (\phi \vec{u}) dV
    
Conditions for applicability:

* :math:`\phi` should be continuous
* :math:`\partial \phi / \partial x_i` should exist
* :math:`\partial \phi / \partial x_i` should be continuous

What is the Hamilton operator in Cartesian coordinates?
=======================================================

Del or Nabla or Hamilton operator:

.. math:: 

    \nabla =  
    \begin{bmatrix}
    \partial / {\partial x_1}  \\
    \partial / {\partial x_2} \\
    \partial / {\partial x_3}
    \end{bmatrix} ^T \cdot 
    \begin{bmatrix}
    \vec{i}  \\
    \vec{j} \\
    \vec{k}
    \end{bmatrix}=
    {\partial \over {\partial x_1}} \vec{i} +
    {\partial \over {\partial x_2}} \vec{j} +
    {\partial \over {\partial x_3}} \vec{k} =
    {\partial \over {\partial x_j}} \vec{e}_j = 
    {\partial \over {\partial x}} \vec{i} + 
    {\partial \over {\partial y}} \vec{j} +
    {\partial \over {\partial z}} \vec{k}
    
What is the Laplacian operator?
===============================

.. math:: 

    \nabla^2 =
    \nabla \cdot \nabla = 
    \begin{bmatrix}
    \partial / {\partial x_1}  \\
    \partial / {\partial x_2} \\
    \partial / {\partial x_3}
    \end{bmatrix} ^T \cdot 
    \begin{bmatrix}
    \partial / {\partial x_1}  \\
    \partial / {\partial x_2} \\
    \partial / {\partial x_3}
    \end{bmatrix}=
    {\partial^2 \over {\partial x_1^2}} +
    {\partial^2 \over {\partial x_2^2}} +
    {\partial^2 \over {\partial x_3^2}} =
    {{\partial^2} \over {\partial x_i^2}} =
    {\partial^2 \over {\partial x^2}} + 
    {\partial^2 \over {\partial y^2}} +
    {\partial^2 \over {\partial z^2}}
    
* Some authors use :math:`\Delta` for the Laplacian

What is the Laplacian of a scalar field?
========================================

.. math:: 

    \nabla^2 u =
    \begin{bmatrix}
    \partial / {\partial x_1}  \\
    \partial / {\partial x_2} \\
    \partial / {\partial x_3}
    \end{bmatrix} ^T \cdot 
    \begin{bmatrix}
    \partial / {\partial x_1}  \\
    \partial / {\partial x_2} \\
    \partial / {\partial x_3}
    \end{bmatrix} u =
    {{\partial^2 u} \over {\partial x_1^2}} +
    {{\partial^2 u} \over {\partial x_2^2}} +
    {{\partial^2 u} \over {\partial x_3^2}} =
    {{\partial^2 u_i} \over {\partial x_j^2}}
    
* :math:`i` = principal direction
* :math:`j` = number of dimensions
* The extend to which the scalr field represents a source (+ve) or sink (-ve)
* In the Navier-Stokes equations, :math:`\nabla^2 \vec{u}` is sometimes used implying :math:`\vec{u}` is a row vector for each principal direction.

What is a rank 0, 1 and 2 tensor?
=================================

* Rank 0 = scalar e.g. :math:`p = p` (components in no direction)
* Rank 1 = vector e.g. :math:`\vec{u} = u_j` (components in one direction)
* Rank 2 = matrix e.g. :math:`\overset{\underset{\mathrm{\rightrightarrows}}{}}{\tau} = \tau_{ij}` (components in two directions)

Imagine a cube:

* :math:`i` = number of directions (columns)
* :math:`j` = principal direction (rows)

.. math::
    
    p = p

.. math::

    u_j = \begin{bmatrix}
    \partial / {\partial u_1}  \\
    \partial / {\partial u_2} \\
    \partial / {\partial u_3}
    \end{bmatrix} ^T

.. math::

    \tau_{ij} = \begin{bmatrix}
    \tau_{11} & \tau_{12} & \tau_{13}  \\
    \tau_{21} & \tau_{22} & \tau_{23} \\
    \tau_{31} & \tau_{32} & \tau_{33}
    \end{bmatrix}
    
What is the tensor product of two vectors?
==========================================

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{T} =
    \vec{u} \otimes \vec{u} = 
    \begin{bmatrix}
    u_1^2 & u_1 u_2 & u_1 u_3 \\
    u_2 u_1 & u_2^2 & u_2 u_3 \\
    u_3 u_1 & u_3 u_2 & u_3^2 
    \end{bmatrix} = 
    u_i u_j

* :math:`i` = row
* :math:`j` = column
    
What is the tensor product of the Hamliton operator and a vector?
=================================================================

.. math::

    \overset{\underset{\mathrm{\rightrightarrows}}{}}{T} =
    \nabla \otimes \vec{u} = 
    \begin{bmatrix}
    {\partial / \partial x} (u) & {\partial / \partial x} (v) & {\partial / \partial x} (w) \\
    {\partial / \partial y} (u) & {\partial / \partial y} (v) & {\partial / \partial y} (w) \\
    {\partial / \partial z} (u) & {\partial / \partial z} (v) & {\partial / \partial z} (w) 
    \end{bmatrix} = 
    \nabla_i u_j
    
* :math:`i` = row
* :math:`j` = column    

What is the divergence of a rank 2 tensor?
==========================================

.. math::

    \vec{u} = 
    \nabla \cdot \overset{\underset{\mathrm{\rightrightarrows}}{}}{T} =
    \begin{bmatrix}
    {\partial / \partial x_1} (T_{11}) + {\partial / \partial x_2} (T_{12}) + {\partial / \partial x_3} (T_{13}) \\
    {\partial / \partial x_1} (T_{21}) + {\partial / \partial x_2} (T_{22}) + {\partial / \partial x_3} (T_{23}) \\
    {\partial / \partial x_1} (T_{31}) + {\partial / \partial x_2} (T_{32}) + {\partial / \partial x_3} (T_{33}) 
    \end{bmatrix} = 
    {\partial \over {\partial x_j}} T_{ij}
    
What is the transpose of the tensor product of the Hamilton operator and a vector?
==================================================================================

.. math::

    (\nabla \otimes \vec{u})^T = \vec{u} \otimes \nabla
    
What is the difference between a symmetric and non-symmetric tensor?
====================================================================

* Symmetric :math:`\longrightarrow` :math:`\vec{a} \otimes \vec{b} = (\vec{a} \otimes \vec{b})^T = \vec{b} \otimes \vec{a}` 
* Non-symmetric :math:`\longrightarrow` :math:`\nabla \otimes \vec{u} = (\nabla \otimes \vec{u})^T \ne \vec{u} \otimes \nabla` 

How do grad, div, curl and the tensor product affect the rank?
==============================================================

.. list-table::
   :header-rows: 1
   :widths: 15 15 50

   * - Symbol
     - Name
     - Meaning
   * - :math:`\nabla_j p`
     - Vector (row)
     - Grad operation increases rank by 1
   * - :math:`\nabla_j u_i`
     - Matrix
     - Grad operation increases rank by 1
   * - :math:`\nabla_j \cdot u_j`
     - Scalar
     - Divergence operator reduces rank by 1
   * - :math:`\nabla_j \cdot T_{ij}`
     - Vector (row)
     - Divergence operator reduces rank by 1
   * - :math:`\nabla \times \vec{u}`
     - Vector
     - Curl operator maintains rank if :math:`\vec{u}` is 3D, or reduces rank by 1 if :math:`\vec{u}` is 2D
   * - :math:`\nabla_i \otimes u_j`
     - Matrix
     - Tensor product increases rank by 1