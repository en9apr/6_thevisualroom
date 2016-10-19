============================================
2D Finite Volume Method: Non-Cartesian Grids
============================================

.. contents::
   :local:

Summary
=======

* The cross product is used to compute the volume
* The midpoint rule is used to compute the flow quantity and the source term
* The FVM with the cell-centred or cell-vertex formulation can be used to compute the fluxes
* If the FVM is applied on a Cartesian grid (every angle is a right angle), the FD formulas are recovered

Non-Cartesian Grid
==================

Non-Cartesian, although shown as Cartesian, but we will not assume we know the volume of the cell ABCD

.. figure:: ../_images/centred_scheme_2.png
   :align: center
   :scale: 70%

For cell ABCD:

.. math:: {\partial \over \partial t} \int_{\Omega_{ij}} U d \Omega + \oint_{ABCD} \mathbf{F} \cdot d \mathbf{S} = \int_{\Omega_{ij}} Q d \Omega

.. math:: {\partial \over \partial t}  U_{ij} d \Omega_{ij} + \sum_{ABCD} (\pm fdy \pm gdx) =  Q_{ij} d \Omega_{ij}

where: :math:`f` and :math:`g` are the Cartesian components of the flux vector :math:`\mathbf{F}`

The sign in the flux appears due to the surface vector, e.g for Side AB:

.. math:: \mathbf{S}_{AB} = \Delta y_{AB} \mathbf{i} - \Delta x_{AB} \mathbf{j} = (y_B - y_A) \mathbf{i} - (x_B - x_A) \mathbf{j}

Evaluation of the Flow Quantity Terms
=====================================

Can use **midpoint rule** which approximates a volume integral by the product of the centre value and the CV:

.. math:: U_{i,j} = \int_{\Omega_{i,j}} U d \Omega = U_{i,j} \Delta \Omega_{i,j}

Evaluation of the Source Terms
==============================

Can use **midpoint rule** which approximates a volume integral by the product of the centre value and the CV:

.. math:: Q_{i,j} = \int_{\Omega_{i,j}} Q d \Omega = Q_{i,j} \Delta \Omega_{i,j}

Evaluation of the Volumes
=========================

* The vector product is the area of a parallelogram

* Take the absolute value to ensure the value is positive (:math:`\Delta x` and :math:`\Delta y` can be negative)

.. math:: \left| \mathbf{a} \times \mathbf{b} \right| = 
          \left| \text{det} \begin{bmatrix} \Delta x_1 & \Delta x_2 \\ \Delta y_1 & \Delta y_2  \end{bmatrix} \right| =
          \left| \Delta x_1 \Delta y_2 - \Delta x_2 \Delta y_1 \right|

* The area we are looking for (the quadrilateral) is half the area of the parallelogram:

.. math:: \Omega_{ABCD} = {1 \over 2} \left| \mathbf{AC} \times \mathbf{BD} \right|

.. math:: \mathbf{AC} = \begin{bmatrix} \Delta x_{AC} \\ \Delta y_{AC} \end{bmatrix}

.. math:: \mathbf{BD} = \begin{bmatrix} \Delta x_{BD} \\ \Delta y_{BD} \end{bmatrix}

.. math:: \Omega_{ABCD} = {1 \over 2} \left|  \Delta x_{AC}  \Delta y_{BD} - \Delta x_{BD} \Delta y_{AC} \right|

.. figure:: ../_images/parallelogram.png
   :align: center
   :scale: 70%

Evaluation of the Fluxes
========================

The evaluation of fluxes determines the difference between a variety of schemes, e.g. cell vertex or cell centered

It depends on the location of the flow variables w.r.t. the mesh and selected scheme

Central Scheme and Cell Centered FVM
------------------------------------

We know the values at the black dots, we don't know the values in red

.. figure:: ../_images/centred_scheme_3.png
   :align: center
   :scale: 70%

* **Option 1:** Average of fluxes perpendicular to face - e.g. **midpoint rule**

.. math:: f_{AB} = f_{i+1/2,j}

Take the average flux:

.. math:: f_{i+1/2,j} =  {1 \over 2} (f_{i,j} + f_{i+1,j})

where :math:`f_{i,j} = f(u_{i,j})`

Flux of the average flow quantity:

.. math:: f_{AB} = f \left( {u_{i,j} + u_{i+1,j} \over 2} \right)

(not the same, because f is generally non-linear)

Can also use Upwind, Linear Interpolation, QUICK and Higher Order Schemes for midpoint rule

* **Option 2:** Average of fluxes parallel to face A and B - e.g. **trapezoid rule**

.. math:: f_{AB} = {1 \over 2} (f_A + f_B)

Evaluate flow quantity in A and B (at the corner A)

.. math:: u_A = {1 \over 4}(u_{i,j} + u_{i+1,j} + u_{i+1,j-1} + u_{i,j-1})

Or average the fluxes:

.. math:: f_A = {1 \over 4}(f_{i,j} + f_{i+1,j} + f_{i+1,j-1} + f_{i,j-1})

(more flux evaluations - could be expensive computationally)

Can also use Simpson's Rule, cubic polynomials for this type of rule

Central scheme and Cell Vertex FVM
----------------------------------

We know the values at the black dots, we don't know the values in red - remember flux is not at a point, it's **through a face** by integration.

.. figure:: ../_images/vertex_scheme_4.png
   :align: center
   :scale: 70%

* **Option 1** - fluxes perpendicular to face:

.. math:: f_{AB} = f \left( {u_{i,j} + u_{i+1,j} \over 2} \right)

* **Option 2** - fluxes parallel to face:

.. math:: f_{AB} = {1 \over 2} (f_A + f_B)

* The last one corresponds to trapezoidal rule for the integral:

.. math:: \int_{AB} f dy = {1 \over 2}(f_A + f_B)(y_B - y_A)

Summing over the sides ABCD:

.. figure:: ../_images/vertex_scheme_6.png
   :align: center
   :scale: 70%

Use trapezium rule: half the sum of the parallel sides times the distance between them:

.. math:: \oint_{ABCD} \mathbf{F} \cdot d \mathbf{S} = 
          {1 \over 2} \left( \begin{bmatrix} {f_A - f_C \\ g_C - g_A} \end{bmatrix} \cdot
                      \begin{bmatrix} {\Delta y_{DB} \\ \Delta x_{DB}} \end{bmatrix} + 
                      \begin{bmatrix} {f_B - f_D \\ g_B - g_D} \end{bmatrix} \cdot
                      \begin{bmatrix} {\Delta y_{AC} \\ \Delta x_{AC}} \end{bmatrix} \right) \\
                      = {1 \over 2} \left( (f_A-f_C) \Delta y_{DB} + (g_C-g_A) \Delta x_{DB} +
                                      (f_B-f_D) \Delta y_{AC} + (g_B-g_D) \Delta x_{AC} \right)



2D Example - FVM Reduces to FD
==============================

This is a generic example - we could use cell centered or cell vertex here (but just showing cell centered):

.. figure:: ../_images/centred_scheme_6.png
   :align: center
   :scale: 70%

For cell ABCD:

.. math:: {\partial \over \partial t} \int_{\Omega_{i,j}} U d \Omega + \oint_{ABCD} \mathbf{F} \cdot d \mathbf{S} = \int_{\Omega_{i,j}} Q d \Omega

.. math:: {\partial \over \partial t} u_{i,j} \Omega_{i,j} +  
           \begin{bmatrix} {f_{AB} \\ g_{AB}} \end{bmatrix} \cdot
           \begin{bmatrix} {\Delta y_{AB} \\ \Delta x_{AB}} \end{bmatrix} + 
           \begin{bmatrix} {f_{BC} \\ g_{BC}} \end{bmatrix} \cdot
           \begin{bmatrix} {\Delta y_{BC} \\ \Delta x_{BC}} \end{bmatrix} +
           \begin{bmatrix} {f_{CD} \\ g_{CD}} \end{bmatrix} \cdot
           \begin{bmatrix} {\Delta y_{CD} \\ \Delta x_{CD}} \end{bmatrix} +
           \begin{bmatrix} {f_{DA} \\ g_{DA}} \end{bmatrix} \cdot
           \begin{bmatrix} {\Delta y_{DA} \\ \Delta x_{DA}} \end{bmatrix} =
           q_{i,j} d \Omega_{i,j}

Volume:

:math:`\Omega_{i,j} = \Delta x \Delta y`


Surfaces:

:math:`\Delta y_{AB}` = :math:`\Delta y`  :math:`\quad`  :math:`\Delta x_{AB}` = :math:`0`  :math:`\quad`  :math:`\Delta y_{BC}` = :math:`0`  :math:`\quad`  :math:`\Delta x_{BC}` = :math:`\Delta x`  :math:`\quad`  :math:`\Delta y_{CD}` = :math:`\Delta y`  :math:`\quad`  :math:`\Delta x_{CD}` = :math:`0`  :math:`\quad`  :math:`\Delta y_{DA}` = :math:`0`  :math:`\quad`  :math:`\Delta x_{DA}` = :math:`\Delta x`

Fluxes:

:math:`f_{AB}` = :math:`f_{i+1/2,j}`  :math:`\quad`  :math:`g_{AB}` = :math:`0`  :math:`\quad`  :math:`f_{BC}` = :math:`0`  :math:`\quad`  :math:`g_{BC}` = :math:`g_{i,j+1/2}`  :math:`\quad`  :math:`f_{CD}` = :math:`-f_{i-1/2,j}`  :math:`\quad`  :math:`g_{CD}` = :math:`0`  :math:`\quad`  :math:`f_{DA}` = :math:`0`  :math:`\quad`  :math:`g_{DA}` = :math:`-g_{i,j-1/2}`

Hence:

.. math:: {\partial \over \partial t} u_{i,j} \Delta x \Delta y + 
          ({f_{i+1/2,j} - f_{i-1/2,j}}) \Delta y + ({g_{i,j+1/2} - g_{i,j-1/2}}) \Delta x = 
          q_{i,j} \Delta x \Delta y

Dividing by the volume we have recovered the finite difference formula:

.. math:: {\partial \over \partial t} u_{i,j} + 
          {({f_{i+1/2,j} - f_{i-1/2,j}}) \over \Delta x} + {({g_{i,j+1/2} - g_{i,j-1/2}}) \over \Delta y} = 
          q_{i,j}

Note: :math:`f_{i,j}` and :math:`g_{i,j}` do not appear, hence **odd-even decoupling** can occur
