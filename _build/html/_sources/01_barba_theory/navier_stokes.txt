===========================================
 Derivation of the Navier-Stokes Equations
===========================================

Assumptions
+++++++++++

* Newtonian fluid: linear relationship between stress and strain, hence viscosity is constant 
* Flow is incompressible: density is constant
* Flow is isothermal: temperature is constant

Normal Stresses
+++++++++++++++

.. math:: \sigma_{xx} = -p + 2 \mu {\partial u \over \partial x}
.. math:: \sigma_{yy} = -p + 2 \mu {\partial v \over \partial y}
.. math:: \sigma_{zz} = -p + 2 \mu {\partial w \over \partial z}

:math:`\nabla \cdot \vec V = 0` for incompressible flow:
````````````````````````````````````````````````````````

.. math:: {1 \over 3} ( \sigma_{xx} + \sigma_{yy} + \sigma_{zz}) = -p

Shear Stresses
++++++++++++++

.. math:: \tau_{xy} = \tau_{yx} = \mu \left ({\partial u \over \partial y} + {\partial v \over \partial x} \right )
.. math:: \tau_{yz} = \tau_{zy} = \mu \left ({\partial v \over \partial z} + {\partial w \over \partial y} \right )
.. math:: \tau_{zx} = \tau_{xz} = \mu \left ({\partial w \over \partial x} + {\partial u \over \partial z} \right )

x-direction Momentum Equation for a Newtonian, incompressible, isothermal fluid
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

From the Momentum Equation:
```````````````````````````

.. math:: \rho g_x + {\partial \sigma_{xx} \over \partial x} + {\partial \tau_{yx} \over \partial y} + {\partial \tau_{zx} \over \partial z} = \rho \left ( {\partial u \over \partial t} + u{\partial u \over \partial x} + v {\partial u \over \partial y} + w {\partial u \over \partial z} \right )

Substitute in normal stresses and shear stresses:
`````````````````````````````````````````````````

.. math:: \rho g_x -{\partial p \over \partial x}+2 \mu {\partial^2 u \over \partial x^2} + \mu \left ({\partial^2 u \over \partial y^2} + {\partial^2 v \over \partial y \partial x} \right ) + \mu \left ({\partial^2 w \over \partial z \partial x} + {\partial^2 u \over \partial z^2} \right ) = \rho \left ( {\partial u \over \partial t} + u{\partial u \over \partial x} + v {\partial u \over \partial y} + w {\partial u \over \partial z} \right )

Rewrite:
````````

.. math:: \rho g_x -{\partial p \over \partial x}+\mu \left ({\partial^2 u \over \partial x^2} + {\partial^2 u \over \partial y^2} + {\partial^2 u \over \partial z^2} \right ) + \mu {\partial \over \partial x} \left ({\partial u \over \partial x} + {\partial v \over \partial y} + {\partial w \over \partial z} \right ) = \rho \left ( {\partial u \over \partial t} + u{\partial u \over \partial x} + v {\partial u \over \partial y} + w {\partial u \over \partial z} \right )

:math:`\nabla \cdot \vec V = 0` for incompressible flow:
````````````````````````````````````````````````````````

.. math:: \rho g_x -{\partial p \over \partial x}+\mu \left ({\partial^2 u \over \partial x^2} + {\partial^2 u \over \partial y^2} + {\partial^2 u \over \partial z^2} \right ) = \rho \left ( {\partial u \over \partial t} + u{\partial u \over \partial x} + v {\partial u \over \partial y} + w {\partial u \over \partial z} \right )

Vector Notation
```````````````

.. math:: \rho \vec g-  \nabla p + \mu \nabla^2 \vec V = \rho \left ( {\partial \vec V \over \partial t} + \vec V (\nabla \cdot \vec V) \right )

Solutions
+++++++++

* 3 Momentum Equations + Continuity = 4 Equations
* Unknowns = u, v, w, p, :math:`\rho = 5` Unknowns
* Need an equation of state - to relate pressure and density
* The Navier-Stokes Equations are time-dependent, non-linear, 2nd order PDEs - very few known solutions (parallel plates, pipe flow, concentric cylinders).

  - The Transient Term is :math:`{\partial \vec V / \partial t}`
  - The Convection Term is :math:`\vec V(\nabla \cdot \vec V)`. This is the non-linear term and is the cause most of the difficulty in solving these equations.
  - The Diffusion Term is :math:`\mu \nabla^2 \vec V`. This is the 2nd order term.
  - The Body Force Term is :math:`\rho \vec g`
  - The Pressure Gradient Term is :math:`\nabla p`

* The only general approach is computational

