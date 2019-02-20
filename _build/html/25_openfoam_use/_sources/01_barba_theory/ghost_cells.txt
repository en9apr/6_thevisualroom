==============================
Ghost Cell Boundary Conditions
==============================

.. contents::
   :local:

Boundary Conditions for Incompressible Navier Stokes
====================================================

In dimensionless form, the Navier-Stokes Equations are:

.. math:: \nabla \cdot \mathbf{u} = 0
   :label: 1

.. math:: {\partial \mathbf{u} \over \partial t} + \mathbf{u} \cdot \nabla \mathbf{u}=
          -\nabla p + {1 \over Re} \nabla^2 \mathbf{u}
   :label: 2

The independent variable is :math:`\mathbf{x} \in \Omega`

**Initial Conditions** are (some function for the initial conditions), with the the velocity divergence free in the domain :math:`\Omega` i.e. :math:`\nabla \cdot \mathbf{u}_0 = 0`

.. math:: \mathbf{u}(\mathbf{x},t_0) = \mathbf{u}_0(\mathbf{x})
   :label: 3

**Boundary Conditions** are at :math:`\mathbf{x} \in \Gamma = \partial \Omega` 

.. math:: \mathbf{u}(\mathbf{x},t) = \mathbf{u}_\Gamma (\mathbf{x},t)
   :label: 4

**Global Mass Conservation** these are at the inflow or outflow boundaries (:math:`\mathbf{n}` is the vector normal to the boundary :math:`\Gamma`)

.. math:: \oint_\Gamma \mathbf{n} \cdot \mathbf{u}_\Gamma d \sigma = 0
   :label: 5

In practice, B.C.s are stipulated separately for the different types:

* inflow
* outflow
* solid surface
* far field

.. figure:: ../_images/boundary_conditions.png
   :scale: 75%
   :align: center

* No slip at body surface :math:`\Rightarrow \mathbf{u}_\Gamma = \mathbf{u}_{body}` (=0 usually)
* Pressure at surface is not known - no equation for p (we have to derive one - e.g. Poisson Equation) - PPE (Pressure Poisson Equation) - the divergence of convection

.. math:: \nabla^2 p = -\nabla \cdot (\mathbf{u} \cdot \nabla \mathbf{u})
  :label: 6

Need to use :eq:`2` and :eq:`6` (:eq:`6` ensures continuity, as the divergence of velocity has been set to zero in it's derivation)

Still don't know BCs for pressure, we have replaced the Continuity Equation (:eq:`1`) with the Poisson Equation (:eq:`6`), so we have increased the order of the equations by one, which means we would need an additional boundary condition, which we don't have.

Fractional Step Method or Projection Method
===========================================

Method used by Rempfer (2006)
-----------------------------

* **Step 1:** Omit pressure term in NS and find intermediate velocity - for explicit formulation (we could use implicit):

.. math:: {{\mathbf{u}^* - \mathbf{u}^n} \over {\Delta t}} = -(\mathbf{u}^n \cdot \nabla) \mathbf{u}^n + {1 \over Re} \nabla^2 \mathbf{u}^n
   :label: 7

* **Step 2:** Apply no slip boundary condition:

.. math:: \left . \mathbf{u}^* \right |_\Gamma = \mathbf{u}_\Gamma (\mathbf{x}, t^{n+1}) \quad \mathbf{x} \in \Gamma
   :label: 8

* **Step 3:** :math:`\mathbf{u}^*` is not necessarily divergence free, so restore divergence free velocity by introducing a "projector" :math:`\tilde{p}`: (this is the projection in the **Pressure Projection Method**)

.. math:: {{ \mathbf{u}^{n+1} - \mathbf{u}^* } \over {\Delta t}} = -\nabla \tilde{p}^{n+1} \quad \mathbf{x} \in \Gamma
   :label: 9

* **Step 4:** Apply Continuity Equation at :math:`n+1` (divergence free condition):

.. math:: \nabla \cdot \mathbf{u}^{n+1} = 0
   :label: 10

* **Step 5:** Apply boundary condition using normal vector:

.. math:: \left . \mathbf{n} \cdot \mathbf{u}^{n+1} \right |_\Gamma = \mathbf{n} \cdot \mathbf{u}_\Gamma (\mathbf{x},t^{n+1}) 
   :label: 11
 
**Observations:** :math:`\tilde{p}` is not identical with pressure - we have taken two steps so the equations contain an error. The pressure is a **pseudo pressure** (numerical pressure) - all it does is correct the velocity so it's divergence free

.. math:: \nabla^2 \tilde{p}^{n+1} ={ {\nabla \cdot \mathbf{u}^*} \over \Delta t}
   :label: 12

Kim and Moin (1984) Say that the pressure, p is given by:

.. math:: p = \tilde{p} + {{\Delta t} \over {2 Re}} \nabla^2 \tilde{p}

Equation for :math:`\tilde{p}`

.. math:: \nabla^2 \tilde{p}^{n+1} = {{\nabla \cdot \mathbf{u}^*} \over {\Delta t}}

By considering BCs :eq:`8`  :eq:`11`  :math:`\Rightarrow \quad \left . \mathbf{n} \cdot \nabla \tilde{p}^{n+1} \right |_\Gamma = 0`

Or:

.. math:: {{\partial \tilde{p}^{n+1}} \over {\partial n}} = 0 
   :label: 13

i.e. homogeneous Neumann Boundary Condition

**Highly controversial issue**

Method by Gresho and Sani (1987)
--------------------------------

Another paper has an alternate method - take the normal of the momentum equation at the boundary to set pressure normal to boundary, everywhere else use Poisson Equation, Continuity Equation and ordinary Momentum Equation, i.e.

.. math:: {\partial u \over \partial t} + \nabla P = \nu \nabla^2 u - u \cdot \nabla u = f

and

.. math:: \nabla \cdot ({\partial u \over \partial t}) = 0 \qquad \text{in} \ \Omega

these imply

.. math:: \nabla^2 P = \nabla \cdot f \qquad \text{in} \ \Omega

also

.. math:: {\partial P \over \partial n} = \mathbf{n} \cdot (f - ({\partial u \over \partial t})) \qquad \text{on} \ \Gamma

Again this gives a Neumann boundary condition - **but Rempfer says that this is not correct**

It looks like one method says that the velocity should be normal to the boundary and the other that pressure should be normal to the boundary.

Method by Kim and Moin
----------------------

Say that pressure boundary conditions are not important because of the use of a staggered grid - **could be naive**

Ghost Cells
===========

Method by Harlow and Welsh 
--------------------------

* Use of "Ghost Cells" - for reflection of field variables near a wall.
* Normal pressure is not a condition set directly at the wall, but implied by the values either side of the wall.
* The same is true for velocity.

Hence, in this scheme there is no difference between the opinions of Rempfer (2006) and Gresho and Sani (1987), although the opinion of Kim and Moin (1984) is clearly slightly misleading.

Boundary Conditions at Rigid Walls
----------------------------------

Rigid walls may be of two types:

* No-slip
* Free-slip (could be a plane of symmetry - or a greased surface)

.. figure:: ../_images/wall.png
   :scale: 200%
   :align: center

* At the wall :math:`v = 0`
* On the fluid side :math:`v \ne 0`
* On the outside :math:`v' = -v` (no slip)
* On the outside :math:`v' = +v` (free slip)

Harlow and Welch use these velocity BCs in order to set the pressure BCs using **"Ghost Cells"**. The velocity at the wall is **not set directly**, but is implied by the pressure BCs.

Free slip wall
~~~~~~~~~~~~~~

**Vertical wall**: from the momentum equation for :math:`u` (with the reversal of all normal velocities and no change in the tangential velocity):

For a wall at :math:`i+{1 \over 2}`, the implied wall velocity :math:`\forall j` is:

.. math:: u_{i+{1 \over 2}} = 0

In the Ghost Cell, the velocities are:

.. math:: u_{i+{3 \over 2}} = -u_{i+{1 \over 2}}  

.. math:: u_{i+{1}} = -u_{i}  

Hence, the pressure BCs to set is:

.. math:: \psi' = \psi \pm g_x \Delta x

:math:`+ \Rightarrow` fluid is on the left of the wall

:math:`- \Rightarrow` fluid is on the right of the wall

Similarly for a **Horizontal wall** at :math:`j + {1 \over 2}`:

.. math:: \psi' = \psi \pm g_y \Delta y

:math:`+ \Rightarrow` fluid is below the wall

:math:`- \Rightarrow` fluid is above the wall

No slip wall
~~~~~~~~~~~~

**Vertical wall**: from the momentum equation for :math:`u` (with the reversal of all tangential velocities and no change in the normal velocity):

For a wall at :math:`i+{1 \over 2}`, the implied wall velocity :math:`\forall j` is:

.. math:: u_{i+{1 \over 2}} = 0

In the Ghost Cell, the velocities are:

.. math:: u_{i+{3 \over 2}} = u_{i+{1 \over 2}}  

.. math:: u_{i+{1}} = u_{i}  

Hence, the pressure BCs to set is:

.. math:: \psi' = \psi \pm g_x \Delta x \pm ({{2 \nu u_{i-{1 \over 2}}} \over {\Delta x}})

:math:`+ \Rightarrow` fluid is on the left of the wall

:math:`- \Rightarrow` fluid is on the right of the wall

Similarly for a **Horizontal wall** at :math:`j + {1 \over 2}`:

.. math::  \psi' = \psi \pm g_y \Delta y \pm ({{2 \nu v_{j-{1 \over 2}}} \over {\Delta y}})

:math:`+ \Rightarrow` fluid is below the wall

:math:`- \Rightarrow` fluid is above the wall





