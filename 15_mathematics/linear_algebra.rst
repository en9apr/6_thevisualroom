=====================================
 Linear Algebra for CFD Applications
=====================================

.. contents::
   :local:

2D Vector Notation
==================

* Divergence of a **vector** gives a scalar (via the scalar product)
* Divergence of a **scalar** doesn't exist
* Gradient of a **vector** gives a vector
* Gradient of a **scalar** gives a scalar

Continuity Equation
-------------------

The 2D Continuity Equation:

----

.. math:: \nabla \cdot \mathbf{u} = 0 
   :label: one

----

|

.. list-table::
   :header-rows: 1
   :widths: 9 9 9 9 9 9

   * - Expression 1
     - Expression 2
     - Component 1a
     - Component 1b
     - Operation
     - Expansion
   * - :math:`\nabla \cdot \mathbf{u}(x,y)`
     - :math:`\text{div } \mathbf{u}`
     - :math:`\begin{bmatrix} \mathbf{i} \  {\partial \over {\partial x}}  \\ \mathbf{j} \ {\partial \over {\partial y}} \end{bmatrix}`
     - :math:`\begin{bmatrix} {u} \\ {v}  \end{bmatrix}`
     - Scalar Product
     - :math:`{{\partial u} \over {\partial x}} + {{\partial v} \over {\partial y}}`

Momentum Equation
-----------------

The 2D Momentum Equation:


----

.. math:: {{\partial \mathbf{u}} \over {\partial t}} + \mathbf{u} \cdot \nabla \mathbf{u} =
          {-{1 \over \rho} \nabla p} + {\nu \nabla^2 \mathbf{u}} 
   :label: two

----

|

.. list-table::
   :header-rows: 1
   :widths: 10 9 10 11 11 9 9 

   * - Expression 1
     - Expression 2
     - Component 1a
     - Component 1b
     - Operation
     - Expansion-x
     - Expansion-y
   * - :math:`{\partial \mathbf{u}(x,y)} \over t`
     - N/A
     - :math:`{\partial \over {\partial t}}`
     - :math:`u` or :math:`v`
     - Time Derivative
     - :math:`{\partial u} \over {\partial t}`
     - :math:`{\partial v} \over {\partial t}`
   * - :math:`\nabla \mathbf{u}(x,y)`
     - :math:`\text{grad } \mathbf{u}`
     - :math:`\begin{bmatrix}  \mathbf{i} \ {\partial \over {\partial x}} \\ \mathbf{j} \ {\partial \over {\partial y}} \end{bmatrix}`
     - :math:`{u}` or :math:`{v}`
     - Vector Gradient
     - :math:`\begin{bmatrix} {{\partial u} \over {\partial x}} \mathbf{i} \\ {{\partial u} \over {\partial y}} \mathbf{j} \end{bmatrix}`
     - :math:`\begin{bmatrix} {{\partial v} \over {\partial x}} \mathbf{i} \\ {{\partial v} \over {\partial y}} \mathbf{j} \end{bmatrix}`
   * - :math:`\mathbf{u} \cdot \nabla \mathbf{u}(x,y)`
     - :math:`\mathbf{u} \cdot \text{grad } \mathbf{u}`
     - :math:`\begin{bmatrix}  {u} \\ {v} \end{bmatrix}`
     - :math:`\begin{bmatrix} {{\partial u} \over {\partial x}} \mathbf{i} \\ {{\partial u} \over {\partial y}} \mathbf{j} \end{bmatrix}` or :math:`\begin{bmatrix} {{\partial v} \over {\partial x}} \mathbf{i} \\ {{\partial v} \over {\partial y}} \mathbf{j} \end{bmatrix}`
     - Scalar Product
     - :math:`u{{\partial u} \over {\partial x}} + v {{\partial u} \over {\partial y}}`
     - :math:`u{{\partial v} \over {\partial x}} + v {{\partial v} \over {\partial y}}`
   * - :math:`\nabla p(x)` or :math:`\nabla p(y)`
     - :math:`\text{grad } p`
     - :math:`{\partial \over {\partial x}}` or :math:`{\partial \over {\partial y}}`
     - :math:`{p}`
     - Scalar Gradient
     - :math:`{{\partial p} \over {\partial x}}`
     - :math:`{{\partial p} \over {\partial y}}`
   * - :math:`\nabla^2 \mathbf{u}(x,y)`
     - :math:`\nabla \cdot \nabla \mathbf{u}(x,y)`
     - :math:`\begin{bmatrix}  \mathbf{i} {\partial \over {\partial x}} \\ \mathbf{j} {\partial \over {\partial y}} \end{bmatrix}`
     - :math:`\begin{bmatrix} {{\partial u} \over {\partial x}} \mathbf{i} \\ {{\partial u} \over {\partial y}} \mathbf{j} \end{bmatrix}` or :math:`\begin{bmatrix} {{\partial v} \over {\partial x}} \mathbf{i} \\ {{\partial v} \over {\partial y}} \mathbf{j} \end{bmatrix}`
     - Laplacian
     - :math:`{{\partial^2 u} \over {\partial x^2}} + {{\partial^2 u} \over {\partial y^2}}`
     - :math:`{{\partial^2 v} \over {\partial x^2}} + {{\partial^2 v} \over {\partial y^2}}`
