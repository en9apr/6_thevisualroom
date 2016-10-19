================================
Finite Difference Formulas in 2D
================================

Extend 1D formulas to 2D
========================

* Just apply the definition of a partial derivative w.r.t. :math:`x` is the variation in :math:`x` holding :math:`y` constant

* Build 2D grid defined by the following:

* For i between :math:`0` and :math:`nx-1`:

.. math::

   x_i = x_0 + i \Delta x

* For j between :math:`0` and :math:`nj-1`:
 
.. math::

   y_j = y_0 + j \Delta y

.. figure:: ../_images/2D.png
   :scale: 75%
   :align: center

* Define :math:`u_{i,j} = u(x_i,y_j)`

Forward Differencing in 2D for 1st derivative
=============================================

* For the point :math:`(i+1,j+1)`, Taylor series in 2D:

.. math:: u_{i+1,j+1} = u_{i,j} +
          \left . \left ( \Delta x {\partial \over \partial x} + \Delta y {\partial \over \partial y} \right )  u \right \vert_{i,j} +
          \left . {1 \over 2} \left ( \Delta x {\partial \over \partial x} + \Delta y {\partial \over \partial y} \right )^2 u \right \vert_{i,j} +
          \left . {1 \over 6} \left ( \Delta x {\partial \over \partial x} + \Delta y {\partial \over \partial y} \right )^3 u \right \vert_{i,j} +
          \cdots + \text{h.o.t}

* 1st order FD in x-direction, i.e.

  - :math:`\Delta y = 0`
  - :math:`j+1=j`
  - Terms of higher order than 1 are zero

.. math::

    \left . {{\partial u} \over {\partial x}} \right \vert_{i,j} = {{u_{i+1,j}-u_{i,j}} \over {\Delta x}} + O(\Delta x)

Central Differencing in 2D for 1st derivative
=============================================

* Take one foward Taylor step, one backward Taylor step 
* **Subtract** the forward and the backward steps
* Re-arrange for the derivative

.. math::

    \left . {{\partial u} \over {\partial x}} \right \vert_{i,j} = {{u_{i+1,j}-u_{i-1,j}} \over {2 \Delta x}} + O(\Delta x^2)

Central Differencing in 2D for 2nd derivative
=============================================

* Take one foward Taylor step, one backward Taylor step 
* **Add** the forward and the backward steps
* Re-arrange for the derivative

.. math::

    \left . {{\partial^2 u} \over {\partial x^2}} \right \vert_{i,j} = {{u_{i+1,j}-2u_{i,j}+u_{i-1,j}} \over {\Delta x^2}} + O(\Delta x^2)
