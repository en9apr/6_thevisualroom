=====================
Second Order Formulas
=====================

.. contents::
   :local:

Question
========

How to define a **second order** finite difference formula for three upwind points? 

Selection of the Mesh Points
============================

Choose: :math:`(i-2)`, :math:`(i-1)` and :math:`i`

Finite Difference Formula with Three Coefficients
=================================================

.. math:: (u_x)_i = {{a u_i + b u_{i-1} + c u_{i-2}} \over {\Delta x}} + O(\Delta x^2)
   :label: 1

Taylor Expansion for i-2 and i-1
================================

.. math::  u_{i-1} = u_i - \Delta x(u_x)_i + 
                    {{(\Delta x)^2} \over 2}(u_{xx})_i -
                    {{(\Delta x)^3} \over 6}(u_{xxx})_i +
                    O(\Delta x^4)
   :label: 2

.. math:: u_{i-2} = u_i - 2 \Delta x(u_x)_i + 
                    {{(2 \Delta x)^2} \over 2}(u_{xx})_i -
                    {{(2 \Delta x)^3} \over 6}(u_{xxx})_i +
                    O(\Delta x^4)
   :label: 3

Substitute :eq:`2` and :eq:`3` into :eq:`1`
===========================================

.. math:: (u_x)_i = {{a u_i + b u_{i-1} + c u_{i-2}} \over {\Delta x}} =
                    {{(a+b+c) u_i - \Delta x(b+2c)(u_x)_i + 
                    {{\Delta x^2} \over 2} (b+4c)(u_{xx})_i} \over {\Delta x}}
   :label: 4

Compare the Coefficients
========================

.. math:: a + b + c = 0 \quad \Rightarrow \quad a = 1.5

.. math:: b + 2c = -1 \quad \Rightarrow \quad -4c + 2c = -1 \quad \Rightarrow \quad c=0.5

.. math:: b + 4c = 0 \quad \Rightarrow \quad b = -4c \quad \Rightarrow \quad b = -2

Update Formula for :math:`u_{x}`
================================

.. math:: (u_{x})_i =  {{3 u_i -4 u_{i-1} + u_{i-2}} \over {2 \Delta x}}

**Notice that the Formula** :math:`a+b+c=0` **applies for all finite difference formulas**

