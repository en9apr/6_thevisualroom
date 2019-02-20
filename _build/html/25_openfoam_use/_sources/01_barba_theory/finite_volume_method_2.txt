===========================
Conservative Discretisation
===========================

.. contents::
   :local:

Definition of a Conservation Law
================================

Recall - what is a conservation law?

"Conservation means the **variation** of a conserved (intensive) flow quantity within a given volume is due to the **net effect** of some **internal sources** and the amount of that quantity which is **crossing the boundary surface**"

Essentially we are looking at the net effect of:

* Sources
* Fluxes

How to Identify Fluxes
----------------------

* Fluxes appear exclusively under the **gradient operator** - the flux is the only spatial derivative that appears in the equation

* This is a way to recognise if flux is conserved - **look at the gradient operator** - can you group the spatial derivatives under a gradient operator, so that the fluxes have **divergence** - i.e. a measure of "outgoingness"?

* In conservative form, we usually have **mass flux**, **momentum flux** and **total energy flux**

How to Identify Sources
-----------------------

* Sources are usually scalars in one direction and don't usually have divergence, e.g. **gravity**, **pressure gradient**

Differential Discretisation
===========================

* The differential form requires the fluxes to be **differentiable** - this is not always satisfied, this may introduce impossible solutions.

* Discontinuities introduce infinite gradient, e.g.

 - shocks in the flow
 - a free surface
 - a slip line behind a trailing edge

Conservative Discretisation
===========================

General Form of a Conservation Law
----------------------------------

.. math:: {{\partial \over \partial t} \int_{\Omega} U d \Omega} + \oint_{S} \mathbf{F} \cdot d\mathbf{S} = \int_{\Omega} Q d \Omega

where:

:math:`\Omega =` volume

:math:`S =` surface

:math:`U =` quantity

:math:`F =` flux vector

:math:`Q =` source

i.e. the quantity :math:`U` only depends on the **surface values** of the fluxes

Conservative Discretisation
---------------------------

:math:`\Omega` divided into :math:`\Omega_1`, :math:`\Omega_2` and :math:`\Omega_3`:

For each subvolume:

.. math:: {{\partial \over \partial t} \int_{\Omega_1} U d \Omega} + \oint_{ABCA} \mathbf{F} \cdot d\mathbf{S} = \int_{\Omega_1} Q d \Omega

.. math:: {{\partial \over \partial t} \int_{\Omega_2} U d \Omega} + \oint_{DEBD} \mathbf{F} \cdot d\mathbf{S} = \int_{\Omega_2} Q d \Omega


.. math:: {{\partial \over \partial t} \int_{\Omega_3} U d \Omega} + \oint_{AEDA} \mathbf{F} \cdot d\mathbf{S} = \int_{\Omega_3} Q d \Omega

Notice that the contribution of the **internal fluxes**, cancel each other out when summated in each direction. Hence the numerical scheme is **conservative** in the way it treats the flux terms.

.. figure:: ../_images/finite_volume_fluxes.png
   :align: center
   :scale: 50%

Example of Conservative Discretisation
--------------------------------------

Consider a 1D conservation law in conservative form:

.. math:: {\partial u \over \partial t} + {\partial f \over \partial x} = q


.. figure:: ../_images/1D_finite_volume.png
   :align: center
   :scale: 70%

Associate a finite "cell" to each mesh point, e.g. cell (i) has two "cell faces" at the mid-points :math:`i-{1/2}` and :math:`i+{1/2}`

Using a CD scheme at point i:

.. math:: {\partial u_i \over \partial t} + {{f_{i+1/2} - f_{i-1/2}} \over \Delta x} = q_i
   :label: 1

Using a CD scheme at point i+1:

.. math:: {\partial u_{i+1} \over \partial t} + {{f_{i+3/2} - f_{i+1/2}} \over \Delta x} = q_{i+1}
   :label: 2

Using a CD scheme at point i-1:

.. math:: {\partial u_{i-1} \over \partial t} + {{f_{i-1/2} - f_{i-3/2}} \over \Delta x} = q_{i-1}
   :label: 3

Add :eq:`1`, :eq:`2` and :eq:`3` and divide by 3 for a consistent discretisation of the conservation law for cell :math:`i-3/2` to :math:`i+3/2`:

.. math:: {\partial \over \partial t} \left( {u_i + u_{i+1} + u_{i-1}} \over 3  \right)  + {{f_{i+3/2} - f_{i-3/2}} \over {3 \Delta x}} = {{q_i + q_{i+1} + q_{i-1}} \over 3}
   :label: 4

A scheme applied directly to cell :math:`i-3/2` to :math:`i+3/2`:

.. math:: {\partial u_i \over \partial t}  + {{f_{i+3/2} - f_{i-3/2}} \over {3 \Delta x}} = {{q_i}}
   :label: 5

The flux of :eq:`4` and :eq:`5` are identical, but the difference is that the flow quantity and source terms are **averaged over three mesh points** in :eq:`4`

Example of Non-Conservative Discretisation
------------------------------------------

Consider a 1D conservation law in conservative form:

.. math:: {\partial u \over \partial t} + a(u){\partial u \over \partial x} = q

where:

:math:`a(u)` is the Jacobian and :math:`a = \partial f / \partial u`

(e.g. Burgers :math:`a=u`, :math:`f=u^2 / 2`)


.. figure:: ../_images/1D_finite_volume_non_conservative.png
   :align: center
   :scale: 70%

Using a CD scheme at point i:

.. math:: {\partial u_i \over \partial t} + a_i{{u_{i+1/2} - u_{i-1/2}} \over \Delta x} = q_i
   :label: 6

Using a CD scheme at point i+1:

.. math:: {\partial u_{i+1} \over \partial t} + a_{i+1}{{u_{i+3/2} - u_{i+1/2}} \over \Delta x} = q_{i+1}
   :label: 7

Using a CD scheme at point i-1:

.. math:: {\partial u_{i-1} \over \partial t} + a_{i-1}{{u_{i-1/2} - u_{i-3/2}} \over \Delta x} = q_{i-1}
   :label: 8

Add :eq:`6`, :eq:`7` and :eq:`8` and divide by 3 for an consistent discretisation of the conservation law for cell :math:`i-3/2` to :math:`i+3/2`:

.. math:: {\partial \over \partial t} \left( {u_i + u_{i+1} + u_{i-1}} \over 3  \right)  + 
     a_i{{u_{i+1/2} - u_{i-1/2}} \over {3 \Delta x}} +
     a_{i+1}{{u_{i+3/2} - u_{i+1/2}} \over {3 \Delta x}} +
     a_{i-1}{{u_{i-1/2} - u_{i-3/2}} \over {3 \Delta x}} =
     {{q_i + q_{i+1} + q_{i-1}} \over 3} +
   :label: 9

Reconstruct what we would want on the LHS, and put the remainder on the RHS:

.. math:: {\partial \over \partial t} \left( {u_i + u_{i+1} + u_{i-1}} \over 3  \right)  + 
     a_i{{u_{i+3/2} - u_{i-3/2}} \over {3 \Delta x}} =
     {{q_i + q_{i+1} + q_{i-1}} \over 3} \\ -
     (a_{i+1}-a_i) {{u_{i+3/2} - u_{i+1/2}} \over {3 \Delta x}} +
     (a_i -a_{i-1}){{u_{i-1/2} - u_{i-3/2}} \over {3 \Delta x}}
   :label: 10

* The extra terms on the RHS of :eq:`10` result from the fact that the flux contributions at the **internal** points of the cell :math:`i-3/2` to :math:`i+3/2` **do not** cancel
* They appear as **additional source terms**
* The discretisation of the non-conservative form of the equation leads to **internal sources**
* The Conservative and Non-Conservative forms of the equations are **mathematically equivalent** for arbitrary non-linear fluxes
* But they are not **numerically equivalent**
* The Taylor expansion of the source terms shows that they are like a 2nd order discretisation proportional to :math:`\Delta x^2(a_x u_x)_x` (same order as the truncation error)
* Numerical experiments show that the non-conservative form is less accurate than the conservative form, especially for discontinuous flows (e.g. shocks) where the gradient is very large, the numerical source terms can become important.



