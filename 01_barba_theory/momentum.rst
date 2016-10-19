=====================================
 Derivation of the Momentum Equation
=====================================

Newton's Second Law - The net force equals the rate of change of momentum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: \vec F = {D {\vec V M} \over D t}

For a system:
+++++++++++++

.. math:: \vec F = {D \over D t} \int_{SYS} \vec V dM

For a control volume (via Reynolds Transport Theorem):
++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
.. math:: \sum \vec F_{CV} = {\partial \over \partial t} \int_{CV} \vec V \rho dV + \int_{CS} \vec V \rho \vec V \cdot \hat n dA 

LHS 1) Body forces - weight for an element :math:`\delta m`:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: \delta \vec F = {D \over D t} \vec V \delta m = \delta m {D \over {Dt}} \vec V = \delta m \cdot \vec a

.. math:: \delta \vec F_b = \delta m \cdot \vec g = \rho \delta x \delta y \delta z \cdot \vec g

LHS 2) Normal force and Tangential force
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../_images/stress_diagram.png

Subscript notation:
  * 1st subscript refers to the direction of the normal vector
  * 2nd subscript refers to the direction of the stress vector

Sign convention:
  * Normal Stress is positive if it's in the same direction as the outward normal vector
  * Shear Stress is positive if it's in the same direction as the coordinate system w.r.t the outward normal vector (i.e. the right hand rule applies - thumb is the direction of the outward normal vector, the two fingers are the shear stresses)

Conservation of Momentum in x-direction:
++++++++++++++++++++++++++++++++++++++++

.. math:: \delta F_{sx} = \left( {\partial \sigma_{xx} \over \partial x} + {\partial \tau_{yx} \over \partial y} + {\partial \tau_{zx} \over \partial z} \right) \delta x \delta y \delta z

Conservation of Momentum in y-direction:
++++++++++++++++++++++++++++++++++++++++

.. math:: \delta F_{sy} = \left( {\partial \tau_{xy}  \over \partial x} + {\partial \sigma_{yy} \over \partial y} + {\partial \tau_{zy} \over \partial z} \right) \delta x \delta y \delta z

Conservation of Momentum in z-direction:
++++++++++++++++++++++++++++++++++++++++

.. math:: \delta F_{sz} = \left( {\partial \tau_{xz} \over \partial x} + {\partial \tau_{yz} \over \partial y} + {\partial \sigma_{zz} \over \partial z} \right) \delta x \delta y \delta z

Equation of motion:
~~~~~~~~~~~~~~~~~~~

.. math:: \rho g_x + {\partial \sigma_{xx} \over \partial x} + {\partial \tau_{yx} \over \partial y} + {\partial \tau_{zx} \over \partial z} = \rho \left ( {\partial u \over \partial t} + u{\partial u \over \partial x} + v {\partial u \over \partial y} + w {\partial u \over \partial z} \right )

	  \rho g_y + {\partial \tau_{xy} \over \partial x} + {\partial \sigma_{yy} \over \partial y} + {\partial \tau_{zy} \over \partial z} = \rho \left ( {\partial v \over \partial t} + u {\partial v \over \partial x} + v {\partial v \over \partial y} + w {\partial v \over \partial z} \right )

	  \rho g_z + {\partial \tau_{xz} \over \partial x} + {\partial \tau_{yz} \over \partial y} + {\partial \sigma_{zz} \over \partial z} = \rho \left ( {\partial w \over \partial t} + u {\partial w \over \partial x} + v {\partial w \over \partial y} + w {\partial w \over \partial z} \right )

Three equations + continuity = Four equations
+++++++++++++++++++++++++++++++++++++++++++++

Unknowns: u, v, w and nine stresses = Twelve unknowns - need more information
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Inviscid flow - no shearing stresses
++++++++++++++++++++++++++++++++++++

.. math:: \sigma_{xx} = \sigma_{yy} = \sigma_{zz}  = -p

Euler's Equation in the x-direction:
++++++++++++++++++++++++++++++++++++

.. math:: \rho g_x - {\partial p \over \partial x} = \rho \left ( {\partial u \over \partial t} + u {\partial u \over \partial x} + v {\partial u \over \partial y} + w {\partial u \over \partial z} \right )

Vector Notation:
++++++++++++++++

.. math:: \rho \vec g- \nabla p = \rho \left ( {\partial \vec V \over \partial t} + \vec V(\nabla \cdot \vec V) \right )
