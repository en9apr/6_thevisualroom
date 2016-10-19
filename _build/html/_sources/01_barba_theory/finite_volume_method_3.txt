=============================================================
2D Finite Volume Method: Spatial and Numerical Discretisation
=============================================================

.. contents::
   :local:

Spatial Discretisation
======================

Types of Spatial Discretisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FVM can handle any type of spatial discretisation

For the same mesh we can have different definitions of the control volumes:

Cell-centred scheme
--------------------

* Control volumes are identical with the grid cells
* Flow variables are located at the centres of grid cells - :math:`q_{i,j}` and :math:`u_{i,j}` can be an average
* Fluxes are located at the volume surfaces (red)

.. figure:: ../_images/centred_scheme.png
   :align: center
   :scale: 70%

Cell-vertex scheme
------------------

* Control volume is 4 cells that share the grid point :math:`i,j`
* Control volumes may overlap

.. figure:: ../_images/vertex_scheme.png
   :align: center
   :scale: 70%

Hexagonal control volume scheme
-------------------------------

.. figure:: ../_images/hexagonal_scheme.png
   :align: center
   :scale: 70%

Trapezoidal control volume scheme
---------------------------------

.. figure:: ../_images/trapezoidal_scheme.png
   :align: center
   :scale: 70%

For a Consistent Scheme
~~~~~~~~~~~~~~~~~~~~~~~

1) The sum of all CVs should cover the entire domain
2) The sub-domains :math:`\Omega_J` are allowed to overlap, as long as each portion of the surface :math:`S_J` have to appear as a part of an **even** number of different surfaces - so that the fluxes cancel out, e.g. 2 surfaces, 4 surfaces etc.
3) Computing fluxes along a cell surface has to be independent of the cell in which they are considered - **ensures that conservation is satisfied**

Numerical Discretisation
========================

Apply integral conservation law to subvolumes (:math:`\Omega_J` = sub-volume, :math:`S_J` = sub-surface)

.. math:: {\partial \over \partial t} \int_{\Omega_J} U d \Omega + 
          \oint_{S_J} \mathbf{F} \cdot d \mathbf{S} = 
          \int_{\Omega_J} Q d \Omega
   :label: 1

Discrete Form, for small volumes (J)

.. math:: {\partial \over \partial t}  (\overline{U_J} \Omega_J) + 
          \sum_{faces} \mathbf{F} \cdot \Delta \mathbf{S_J} = 
          \overline{Q_J} \Omega_J
   :label: 2

Where:

:math:`\overline{U_J}` and :math:`\overline{Q_J}` are averaged values over the cell

If integrating between time level :math:`n` and :math:`n+1`, then **using 1st order Euler** (we could use Runge-Kutta etc)

.. math:: (\overline{U_J} \Omega_J)^{n+1} = (\overline{U_J} \Omega_J)^{n}
          - \Delta t \sum_{faces} \mathbf{F^*} \cdot \Delta \mathbf{S_J} + 
          \Delta t \overline{Q_J^*} \Omega_J
   :label: 3

What are the Cell Averaged Values and the Star Values?

Cell-Averaged Values

.. math:: \overline{U_J^n} = \left. {1 \over \Omega_J} \int_{\Omega_J} U d \Omega_J \right|^n
   :label: 4

.. math:: \overline{Q_J} = \left. {1 \over \Omega_J} \int_{\Omega_J} Q d \Omega_J \right|^n
   :label: 5

Star Values

.. math:: \mathbf{F^*} \cdot \Delta \mathbf{S_J} = {1 \over \Delta t} \int_n^{n+1} \mathbf{F} \cdot \mathbf{\Delta S_J} dt
  :label: 6

.. math:: \overline{Q_J^*} = {1 \over \Delta t} \int_n^{n+1} \overline{Q_J} dt
  :label: 7

The reason we omitted the time index in :eq:`3` for the balance of fluxes and the sources was to indicate the choice:

* If :math:`n` was chosen, it would have been an explicit scheme
* If :math:`n+1` was chosen it would have been an implicit scheme

A scheme is identified by the way the numerical flux :math:`\mathbf{F^*}` approximates the time averaged physical flux across each cell face.

Leaving open the choice of time integrator:

.. math:: {d \over dt} (\overline{U_J} \Omega_J) = -
          \sum_{faces} \mathbf{F^*} \cdot \Delta \mathbf{S_J} + 
          \overline{Q_J^*} \Omega_J = -R_J

RHS defines the "residual" :math:`R_J`
