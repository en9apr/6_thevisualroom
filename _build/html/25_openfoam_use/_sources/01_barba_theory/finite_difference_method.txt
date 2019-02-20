==============================
 The Finite Difference Method
==============================

Origin and Concept
==================

* Invented by Euler in 1768 for one dimension, extended by Runge in 1908 to two dimensions
* Concept is to approximate derivatives using Taylor Expansions

Define numerical grid
=====================
 
.. image:: ../_images/grid.png
   :scale: 75%
   :align: center

* Two families of grid lines
* Grid lines of the same family do not intersect
* Grid lines of different families intersect only once
* Node (i,j) is in 2D - an unknown field variable which depends on neighbouring nodes providing one algebraic equation

Define Derivative
=================

Mathematical Interpretation
---------------------------

.. math:: \left .  {\partial u \over \partial x} \right \vert_i = \lim_{\Delta x \rightarrow 0} {u(x_i + \Delta x) - u(x_i) \over \Delta x}

Geometric Interpretation
------------------------

* Slope of the tangent to the curve with three approximation to the exact solution: Backward, Forward and Central Difference.

.. image:: ../_images/derivatives.png
   :align: center
   :scale: 90%

* Backward difference

.. math:: \left .  {\partial u \over \partial x} \right \vert_i \approx {{u_i - u_{i-1}} \over \Delta x}

* Forward difference

.. math:: \left .  {\partial u \over \partial x} \right \vert_i \approx {{u_{i+1} - u_i} \over \Delta x}

* Central difference

.. math:: \left .  {\partial u \over \partial x} \right \vert_i \approx {{u_{i+1} - u_{i-1}} \over 2 \Delta x}

Error
-----

* Some approximations are better than others
* Quality of approximation improves as :math:`\Delta x` is made smaller

Taylor Series Expansion - Order of the approximations
-----------------------------------------------------

.. math:: u(x) = u(x_i)+(x-x_i) \left . {\partial u \over \partial x} \right \vert_i + {(x - x_i)^2 \over 2!} \left . {\partial^2 u \over \partial x^2} \right \vert_i + \cdots + {(x - x_i)^n \over n!} \left . {\partial^n u \over \partial x^n} \right \vert_i

* Forward differencing: :math:`x = x_{i+1}`
* Backward differencing: :math:`x = x_{i-1}`

Forward Differencing
~~~~~~~~~~~~~~~~~~~~

* We need to obtain the derivative :math:`\left . {\partial u \over \partial x} \right \vert_i`

.. math:: \left . {\partial u \over \partial x} \right \vert_i = {(u_{i+1} - u_i) \over (x_{i+1} - x_i)} -  {(x_{i+1} - x_i) \over 2!} \left . {\partial^2 u \over \partial x^2} \right \vert_i - \cdots - {(x_{i+1} - x_i)^{n-1} \over n!} \left . {\partial^n u \over \partial x^n} \right \vert_i


* If :math:`x_{i+1} - x_i` is small, then:

.. math:: \left . {\partial u \over \partial x} \right \vert_i = {(u_{i+1} - u_i) \over \Delta x} - O(\Delta x)

* There is a possibility that the derivative :math:`\left . {\partial^2 u \over \partial x^2} \right \vert_i` is large, but we assume that the function is well-behaved.

* Forward differencing approximation neglected terms of :math:`O(\Delta x)` :math:`\rightarrow` TRUNCATION ERROR

* As :math:`\Delta  x \rightarrow 0`  :math:`\Rightarrow` FD converges!

Central Differencing
~~~~~~~~~~~~~~~~~~~~

* For Central Differencing, the error is :math:`O(\Delta x^2)`

