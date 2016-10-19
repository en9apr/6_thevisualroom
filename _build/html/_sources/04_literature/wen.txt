=================================
Wen (2012) "Wet/Dry Areas Method"
=================================

Wen, X. "Wet-dry areas method for interfacial (free surface) flow", International Journal for Numerical Methods in Fluids, vol 71, pp 316 - 338, 2012

.. contents::
   :local:

Where is this study placed within the wave modelling community?
===============================================================

Non-linear Breaking Waves in Shallow Water
------------------------------------------

* This paper has an important application to non-linear breaking waves in shallow water. The approach taken in this category follows the fully non-linear equations for a deterministic time-resolving model for one realisation of wave frequency, similar to the approach of Ravel et al (2009). However, unlike Ravel, the numerical method used reformulates the convective mass fluxes in order to allow highly non-linear behaviour at the interface. This allows breaking waves to be considered, which was not considered by Ravel et al (2009). 

* It is known in the wave modelling community that as the water depth reduces, wave breaking becomes an important phenomenon. However, there is also the simultaneous effect of bottom friction, and it is difficult to separate the effect of bottom friction from the effect of wave breaking in shallow water. Although this study focusses mainly on wave breaking, it follows the approach outlined by Cavaleri et al (2007) in which wave breaking and bottom friction could be studied together if needed.

* The time resolving approach may be better suited to non-linear breaking in shallow water compared with stochastic methods because of the highly non-linear and dissipative effects that are present. However, there is an inherent computational penalty of the time-resolving approach as multiple realisations may be needed.

Numerics
--------

* This study is also in the class of numerics, attempting to minimise the numerical errors created due to the description of the flow field.
* The approach has the potential to introduce less numerical error compared to the arithmetic mean approach for the density term in the convective mass fluxes of say Zwart, Burns and Galpin (2007).
* The approach is 2D, such that the wet-dry areas are really wet-dry lengths, but most numerical simulations of this type are 2D.

Summary of Findings
===================

New numerical method for 2D free surface flows using:

* Standard staggered grid
* **Conservative integral form** of Navier-Stokes equations
* 2 continuity equations: 1) mass conservation, 2) volume conservation
* Convective mass flux computed from wet/dry area exposed

Main benefit compared to standard approaches:

* The mean density approach can be shown to produce 450 times the convective mass flux compared with the wet-dry areas method.
* The mean viscosity approach produces 10 times the diffusive mass flux, due to the smaller change in viscosity compared to density across the interface (so this was not considered).

Comparison between this new method shows good agreement with:

* Analytical solution
* Experimental images
* Experimental velocities
      
Literature Review
=================

Problems
--------

1) Instantaneous density jump across a free surface - use conservative integral form of Navier-Stokes because density in convection term stays inside integrand
2) Computing accurate velocity field, hence shear stress for wind blowing over waves - use finite volume method to ensure conservation of mass and momentum

Purposes
--------

**Purpose 1: Explore conservative integral form for free surface flows**

* The conservative integral form of the Navier-Stokes equations has been used in the **cut-cell method** - to cut solid boundaries from the Cartesian grid to include only fluid in the calculation
* High resolution Godunov scheme was used for free surface flows to calculate density and velocity
* Cut cell method has been used to cut air out of the solution

**Purpose 2: Accurately calculate mass fluxes using standard staggered grid**

* Mass fluxes passing through the u and v control volumes are the most important
* For co-located grids, mass fluxes on the f control volumes have been for the mass fluxes on the u and v control volumes before
* For staggered grids, more difficult to get u and v fluxes from f flux - can use an f-grid that is twice as fine as the u and v grid - but this is not standard

**Other features**

* Second order scheme for time derivatives
* High-resolution scheme for interpolation of velocity onto the surface of the control volume
* Second order scheme for viscous terms
* SIMPLEC-PISO scheme for pressure-velocity coupling
* Explicit high resolution compressive interface capturing scheme for arbitrary meshes (CICSAM) for predicting the volume fraction of the fluid

Mathematical Formulation
========================

Governing Equations
-------------------

* Volume Conservation:

.. math:: \int_{S}( \mathbf n \cdot \mathbf v ) dS = 0
   :label: volume

* Mass Conservation:

.. math::  \int_V {\partial \rho \over \partial t} dV +
           \int_S (\rho \mathbf n \cdot \mathbf v) dS = 0
   :label: mass

* Momentum Conservation for :math:`u`:

.. math:: \int_V {\partial (\rho u) \over \partial t} dV +
          \int_S {(\rho \mathbf n \cdot \mathbf v)u} dS = 
          -\int_S {n_x p} dS + 
          \int_S {\mu {\partial u \over \partial n}} dS
   :label: momentum_u

* Momentum Conservation for :math:`v`:

.. math:: \int_V {\partial (\rho v) \over \partial t} dV +
          \int_S {(\rho \mathbf n \cdot \mathbf v)v} dS = 
          -\int_S {n_y p} dS + 
          \int_S {\mu {\partial v \over \partial n}} dS - mg
   :label: momentum_v

* Volume Fraction Conservation:

.. math:: \int_V {\partial f \over \partial t} dV +
           \int_S (\mathbf n \cdot \mathbf v)f dS = 0
   :label: volume_fraction

where:

:math:`\mathbf n = n_x \mathbf i + n_y \mathbf j = \begin{bmatrix} 1 & 0 \end{bmatrix} \text{for east surface}, \begin{bmatrix} 0 & 1 \end{bmatrix} \text{for north surface}, \begin{bmatrix} -1 & 0 \end{bmatrix} \text{for west surface}, \begin{bmatrix} 0 & -1 \end{bmatrix} \text{for south surface}`

:math:`\mathbf v = u \mathbf i + v \mathbf j` = :math:`\begin{bmatrix} u \\ v \end{bmatrix}`

:math:`n_x = 1 \text{ for east surface}, -1 \text{ for west surface}`, :math:`n_y = 1 \text{ for north surface}, -1 \text{ for south surface}`

:math:`{\partial \over \partial n} = {\partial \over \partial x} \text{ for east surface}, -{\partial \over \partial x} \text{ for west surface}`,
:math:`{\partial \over \partial y} \text{ for north surface}, -{\partial \over \partial y} \text{ for south surface}`


Unknowns:

* Velocities :math:`u` and :math:`v`
* Pressure :math:`p`
* Volume fraction :math:`f`

Discretisation of Governing Equations
-------------------------------------

**Aim: Obtain an algebraic formulation for velocities at centre of control volume** :math:`u_P` **and** :math:`v_P`

.. figure:: ../_images/staggered_location.png
   :scale: 75%
   :align: center

**LHS:** Staggered Location for :math:`u, v, p` and :math:`f` **RHS:** :math:`u` control volume

**Discrete Mass Conservation Equation**

Mass Conservation Equation :eq:`mass` - all values are at time n+1

.. math::  {\partial m \over \partial t} +
           \int_{A_e} {\rho u} dS -
           \int_{A_w} {\rho u} dS +
           \int_{A_n} {\rho v} dS -
           \int_{A_s} {\rho v} dS = 0 \\
           \Rightarrow {\partial m \over \partial t} +
           (\rho u A)_e -
           (\rho u A)_w +
           (\rho v A)_n -
           (\rho v A)_s = 0 \\
           \Rightarrow {\partial m \over \partial t} +
           F_e - F_w + F_n - F_s = 0
   :label: mass_int

**Discrete Momentum Conservation Equation for** :math:`u`

Momentum Conservation Equation :eq:`momentum_u` - all values are at time n+1

.. math::  {\partial (u_P m) \over \partial t} +
           \int_{A_e} {(\rho u)u} dS -
           \int_{A_w} {(\rho u)u} dS +
           \int_{A_n} {(\rho v)u} dS -
           \int_{A_s} {(\rho v)u} dS = \\
           \int_{A_w} {p} dS -
           \int_{A_e} {p} dS +
           \int_{A_e} {\mu{\partial u \over \partial x}} dS -
           \int_{A_w} {\mu{\partial u \over \partial x}} dS +
           \int_{A_n} {\mu{\partial u \over \partial y}} dS -
           \int_{A_s} {\mu{\partial u \over \partial y}} dS \\
           \Rightarrow  {\partial (u_P m) \over \partial t} +
           (\rho u A u)_e -
           (\rho u A u)_w +
           (\rho v A u)_n -
           (\rho v A u)_s = \\
           (p A)_w -
           (p A)_e +
           \left ( {\mu A {\partial u \over \partial x}} \right )_e -
           \left ( {\mu A {\partial u \over \partial x}} \right )_w +
           \left ( {\mu A {\partial u \over \partial y}} \right )_n -
           \left ( {\mu A {\partial u \over \partial y}} \right )_s \\
           \Rightarrow  {\partial (u_P m) \over \partial t} +
           F_e u_e - F_w u_w + F_n u_n - F_s u_s = \\
           (p A)_w - (p A)_e +
           D_e(u_E - u_P) - D_w(u_P - u_W) + D_n(u_N - u_P) - D_s(u_P - u_S)
   :label: mom_int

Apply Numerical Schemes
-----------------------

Transient Terms: Second Order Backward Euler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::  {\partial m \over \partial t} = 
           {m_{n+1}-m^{n} \over  \Delta t} +
           {m_{n+1}-2m^{n}+m^{n-1} \over 2 \Delta t} 
   :label: trans_1

.. math::  {\partial (u_Pm) \over \partial t} = 
           {(mu_P)^{n+1}-(mu_P)^{n} \over  \Delta t} +
           {(mu_P)^{n+1}-2(mu_P)^{n}+(mu_P)^{n-1} \over 2 \Delta t}
   :label: trans_2

**Simplification to 1st Order (to illutrate numerical scheme):** all values are at n+1 unless stated otherwise (including the derivative)

.. math::  {\partial m \over \partial t} = 
           {m-m^n \over  \Delta t} 
   :label: trans_3

.. math::  {\partial (u_Pm) \over \partial t} = 
           {(mu_P)-(mu_P)^n \over  \Delta t} 
   :label: trans_4

Velocity Terms at Faces: Blended Scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::  u_e = u_P + \beta {{u_P-u_W} \over 2} 
   :label: face_velocity

where for flow in +ve x-direction: :math:`\beta = 0 \text{ (1st order upwind)}` :math:`\beta = {{u_E-u_P} \over {u_P-u_W}} \text{ (central)}` and :math:`\beta = 1 \text{ (2nd order upwind)}`

For flow in other directions:

================ ============== ============= ========================== ============================== ============================== 
Face velocity    +ve upwind     -ve upwind    Central                    +ve 2nd order upwind           -ve 2nd order upwind
================ ============== ============= ========================== ============================== ==============================
:math:`u_e`      :math:`u_P`    :math:`u_E`   :math:`{u_P+u_E} \over 2`  :math:`{3u_P-u_W} \over 2`     :math:`{3u_E-u_{EE}} \over 2`
:math:`u_w`      :math:`u_W`    :math:`u_P`   :math:`{u_W+u_P} \over 2`  :math:`{3u_W-u_{WW}} \over 2`  :math:`{3u_P-u_E} \over 2`
:math:`u_n`      :math:`u_P`    :math:`u_N`   :math:`{u_P+u_N} \over 2`  :math:`{3u_P-u_S} \over 2`     :math:`{3u_N-u_{NN}} \over 2`
:math:`u_s`      :math:`u_S`    :math:`u_P`   :math:`{u_S+u_P} \over 2`  :math:`{3u_S-u_{SS}} \over 2`  :math:`{3u_P-u_N} \over 2`
================ ============== ============= ========================== ============================== ============================== 

.. figure:: ../_images/discretisation_schemes.png
   :scale: 75%
   :align: center

**Simplification to 1st order upwind scheme (to illustrate numerical scheme)**

:math:`F_e > 0 \Rightarrow u_e = u_P`

:math:`F_e < 0 \Rightarrow u_e = u_E`

.. math::  F_e u_e = u_P \begin{bmatrix}F_e, & 0 \end{bmatrix}  - u_E\begin{bmatrix}-F_e, & 0 \end{bmatrix}
   :label: face_velocity_e

where  :math:`\begin{bmatrix} F_e, & 0 \end{bmatrix}` means the maximum of :math:`F_e` and 0

Similarly:

.. math:: F_w u_w = u_W \begin{bmatrix}F_w, & 0 \end{bmatrix}  - u_P\begin{bmatrix}-F_w, & 0 \end{bmatrix}
   :label: face_velocity_w

.. math:: F_n u_n = u_P \begin{bmatrix}F_n, & 0 \end{bmatrix}  - u_N\begin{bmatrix}-F_n, & 0 \end{bmatrix}
   :label: face_velocity_n

.. math:: F_s u_s = u_S \begin{bmatrix}F_s, & 0 \end{bmatrix}  - u_P\begin{bmatrix}-F_s, & 0 \end{bmatrix}
   :label: face_velocity_s

Apply the First Order Schemes to the Momentum and Continuity Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Rule: With no source terms, and in steady flow, when all the values of** :math:`u` **are the same at the edges, the central value of** :math:`u` **must be equal to them. This means we must take** :math:`u_P` **times the continuity equation away from the momentum equation** 

Momentum Equation :eq:`mom_int` - :math:`u_P` times Continuity Equation :eq:`mass_int`

.. math:: {m_P u_P-m_P^n u_P^n \over  \Delta t} - {m_P u_P -m_P^n u_P \over  \Delta t}
          + u_P \begin{bmatrix}F_e, & 0 \end{bmatrix}  - u_E\begin{bmatrix}-F_e, & 0 \end{bmatrix} - F_e u_P 
          - u_W \begin{bmatrix}F_w, & 0 \end{bmatrix}  - u_P\begin{bmatrix}-F_w, & 0 \end{bmatrix} + F_w u_P \\
          + u_P \begin{bmatrix}F_n, & 0 \end{bmatrix}  - u_N\begin{bmatrix}-F_n, & 0 \end{bmatrix} - F_n u_P 
          - u_S \begin{bmatrix}F_s, & 0 \end{bmatrix}  - u_P\begin{bmatrix}-F_s, & 0 \end{bmatrix} + F_s u_P = \\
           p_w A_w - p_e A_e + D_e(u_E - u_P) - D_w(u_P - u_W) + D_n(u_N - u_P) - D_s(u_P - u_S)

Apply the identities:

:math:`\begin{bmatrix}F_e, & 0 \end{bmatrix} - F_e = \begin{bmatrix}-F_e, & 0 \end{bmatrix}`

:math:`\begin{bmatrix}-F_w, & 0 \end{bmatrix} + F_w = \begin{bmatrix}F_w, & 0 \end{bmatrix}`

:math:`\begin{bmatrix}F_n, & 0 \end{bmatrix} - F_n = \begin{bmatrix}-F_n, & 0 \end{bmatrix}`

:math:`\begin{bmatrix}-F_s, & 0 \end{bmatrix} + F_s = \begin{bmatrix}F_s, & 0 \end{bmatrix}`

.. math:: \left ({m_P^n \over  \Delta t}
          +  D_e + \begin{bmatrix}-F_e, & 0 \end{bmatrix}
          + D_w + \begin{bmatrix}F_w, & 0 \end{bmatrix} 
          + D_n + \begin{bmatrix}-F_n, & 0 \end{bmatrix} 
          + D_s + \begin{bmatrix}F_s, & 0 \end{bmatrix} \right ) u_P = \\
          {m_P^n \over  \Delta t}u_P^n
          + \left (D_e + \begin{bmatrix}-F_e, & 0 \end{bmatrix} \right ) u_E
          + \left (D_w + \begin{bmatrix}F_w, & 0 \end{bmatrix} \right ) u_W
          + \left (D_n + \begin{bmatrix}-F_n, & 0 \end{bmatrix} \right ) u_N
          + \left (D_s + \begin{bmatrix}F_s, & 0 \end{bmatrix} \right ) u_S 
          + p_w A_w - p_e A_e
      

Collect like terms: higher order terms (h.o.t) would be collected on the R.H.S. (all values at n+1 unless otherwise stated):

.. math:: a_P u_P =  {m_P^n \over  \Delta t} u_P^n + a_E u_E + a_W u_W + a_N u_N + a_S u_S + p_w A_w - p_e A_e + \text { h.o.t}
   :label: u_momentum_equation

Similarly for the :math:`v` momentum equation and continuity equation (all values at n+1 unless otherwise stated):

.. math:: a_P v_P =  {m_P^n \over  \Delta t} v_P^n + a_E v_E + a_W v_W + a_N v_N + a_S v_S + p_s A_s - p_n A_n - mg + \text { h.o.t}
   :label: v_momentum_equation

**Note:** The Rule we needed in steady state with no source terms holds, since :math:`a_P = a_E + a_W + a_N + a_S`

Such that in this situation :math:`u_P = u_E = u_W = u_N = u_S` 

**The wet-dry areas method must now DESCRIBE the following for u-control volume and v-control volume:**

* **Convective fluxes** :math:`F_e, F_w, F_n, F_s` (using new method for wet-dry lengths at each face)
* **Diffusive fluxes**  :math:`D_e, D_w, D_n, D_s` (using conventional harmonic mean viscosity approach)
* **Mass of control volume** :math:`m_P` (using conventional arithmetic mean density approach)

The Wet-Dry Areas Method
========================

Terms integrated over the surface:

* Volume flux :math:`m^3s^{-1}m^{-2}` (in the Volume Conservation Equation :eq:`volume` and Volume Fraction Conservation Equation :eq:`volume_fraction`)
* Mass flux :math:`kgs^{-1}m^{-2}` (in the Mass Conservation Equation :eq:`mass`)
* Convective acceleration (in the u Momentum Equation :eq:`momentum_u` and the v Momentum Equation :eq:`momentum_v`)
* Diffusion (in the u Momentum Equation :eq:`momentum_u` and the v Momentum Equation :eq:`momentum_v`)

Terms integrated over the volume:

* Unsteady acceleration (in the Mass Conservation Equation :eq:`mass`, u Momentum Equation :eq:`momentum_u`, v Momentum Equation :eq:`momentum_v` and Volume Fraction Conservation Equation :eq:`volume_fraction`)
* Body forces (in the v Momentum Equation :eq:`momentum_v`)

u-control volume
----------------

**Determine the Flux Terms by considering the Convective Term for the North Face:**

.. math:: \int_{A_n} {(\rho v)u} dS = v_n u_n \int_{L_w+L_a} {\rho dS} =
   v_n u_n \left ( \int_{L_w} {\rho_w dS} + \int_{L_a} {\rho_a dS} \right ) =
   v_n u_n ( \rho_w L_w + \rho_a L_a ) =
   F_n u_n

.. math:: F_n = v_n ( \rho_w L_w + \rho_a L_a )_n

Similarly:

.. math:: F_e = u_e ( \rho_w L_w + \rho_a L_a )_e

.. math:: F_w = u_w ( \rho_w L_w + \rho_a L_a )_w

.. math:: F_s = v_s ( \rho_w L_w + \rho_a L_a )_s


.. figure:: ../_images/wet_dry.png
   :scale: 75%
   :align: center

**What value do the terms** :math:`L_a` **and** :math:`L_w` **have?**

=============== ======================= =======================
State           :math:`L_a`             :math:`L_w`
=============== ======================= =======================
100% air        :math:`\Delta x`        0
100% water      0                       :math:`\Delta x`
air-water mix   :math:`\Delta x - L_w`   :math:`\Delta x - L_a`
=============== ======================= =======================

**Outward Normal Vector**

Negative sign indicates that the change in volume fraction in the positive x and y directions is negative

.. math:: \mathbf n = {{- \nabla \cdot f} \over {\left\vert \nabla \cdot f \right\vert}} = 
  -{{1} \over {\left\vert \nabla \cdot f \right\vert}}
   {\left ( {{\partial f} \over {\partial x}} \mathbf i + {{\partial f} \over {\partial y}} \mathbf j \right )}

* Case 1: :math:`-{{\partial f} \over {\partial y}}>0` and :math:`\left\vert{{\partial f} \over {\partial y}}\right\vert>\left\vert{{\partial f} \over {\partial x}}\right\vert` 

* Case 2: :math:`-{{\partial f} \over {\partial y}}<0` and :math:`\left\vert{{\partial f} \over {\partial y}}\right\vert>\left\vert{{\partial f} \over {\partial x}}\right\vert`

* Case 3: :math:`-{{\partial f} \over {\partial x}}>0` and :math:`\left\vert{{\partial f} \over {\partial x}}\right\vert>\left\vert{{\partial f} \over {\partial y}}\right\vert`

* Case 4: :math:`-{{\partial f} \over {\partial x}}<0` and :math:`\left\vert{{\partial f} \over {\partial x}}\right\vert>\left\vert{{\partial f} \over {\partial y}}\right\vert`

.. figure:: ../_images/normals.png
   :scale: 75%
   :align: center

Now we only need to calculate the wet dry areas case 1, and then rotate it's orientation for the other cases.


:math:`L_w` and :math:`L_a` where :math:`tan(\alpha) < \Delta y / \Delta x` on north face
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Consider Case 1**

Compute for each u and v control volume:

* Wet area
* Dry area
* Viscosity on each face
* Total mass

Assumptions:

* Angle between unit vector :math:`\mathbf n` and positive y axis is in the range -45 degrees to +45 degrees
* Bottom to top, :math:`f` and :math:`\mu` are discontinuous
* Left to right, the :math:`f` is a continous function and conventional interpolation or averaging can be applied
* Air-water interface is a straight line in control volume

.. figure:: ../_images/similar_triangles.png
   :scale: 75%
   :align: center

.. math:: tan(\alpha) = \left \vert {{(f_{i+1,j+1}+f_{i+1,j}-f_{i,j+1}-f_{i,j})\Delta y} \over {(f_{i,j+1}+f_{i+1,j+1}-f_{i,j}-f_{i+1,j})\Delta x}} \right \vert
  :label: alpha

A table might also clarify:

======================== ================ ================== ================== ==================== =================================== ================
Case                     :math:`f_{i,j}`  :math:`f_{i+1,j}`  :math:`f_{i,j+1}`  :math:`f_{i+1,j+1}`  :math:`tan(\alpha)`                 :math:`\alpha`      
======================== ================ ================== ================== ==================== =================================== ================
1                        1                0.5                0.5                0                    :math:`\left \vert 1  \right \vert` 45
2                        0                0.5                0.5                1                    :math:`\left \vert 1  \right \vert` 45
1 flipped horizontally   0.5              1                  0                  0.5                  :math:`\left \vert -1 \right \vert` 45
2 flipped horizontally   0.5              0                  1                  0.5                  :math:`\left \vert -1 \right \vert` 45
======================== ================ ================== ================== ==================== =================================== ================

The volume fraction is continuous from left to right, so the standard average process is applied:

.. math:: f_a = {{f_{i+1,j}+f_{i,j}} \over 2}
   :label: FA

.. math:: f_b = {{f_{i+1,j+1}+f_{i,j+1}} \over 2}
   :label: FB

.. figure:: ../_images/similar_triangles_2.png
   :scale: 75%
   :align: center

:math:`f` is defined as:

.. math::

   f = {{\text{area of water}} \over {\text{area of cell}}}

:math:`1-f` is defined as:

.. math::

   1-f = {{\text{area of air}} \over {\text{area of cell}}}


For similar triangles:

.. math::

   \left ( {{\text{length of water}} \over {\text{length of air}}} \right )^2 = 
   {{\text{area of water}} \over {\text{area of air}}}

Hence:

.. math::
  \left ( {L_w} \over {L_a} \right )^2 = 
   {{f_b \Delta x \Delta y_{j+1}} \over {(1-f_a) \Delta x \Delta y_{j}}}

**Solution for** :math:`L_w`

.. math::  L_w = {({f_b \Delta y_{j+1})^{0.5}} \over {{(f_b \Delta y_{j+1})^{0.5}} + ((1-f_a) \Delta y_{j})^{0.5} }}
   :label: LW

**Solution for** :math:`L_a`

.. math:: L_a = \Delta x - L_w
   :label: LA

From Equation :eq:`LW` and :eq:`LA`

* If :math:`f_b = 0` and :math:`f_a = 1` then :math:`L_w = 0` and :math:`L_a = \Delta x`
* If :math:`f_b = 1` and :math:`f_a = 0` then :math:`L_w = \Delta x` and :math:`L_a = 0`

Advantages of Equation :eq:`LW`:

* We don't need to know that the interface is between P and N for Equation :eq:`LW`
* We don't need the angle :math:`\alpha`

:math:`L_w` and :math:`L_a` where :math:`tan(\alpha) > \Delta y / \Delta x` on north face
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Case A solution for** :math:`L_w` **(** :math:`L_a` **from equation** :eq:`LA` **)**

.. math:: f_a V_a + f_b V_b \leqslant {{(\Delta y_j)^2} \over {2 tan \alpha}} \\

.. math:: L_w = 0
   :label: LWA

.. math:: f_b = 0

.. figure:: ../_images/limits_a.png
   :scale: 75%
   :align: center

**Case B solution for** :math:`L_w` **(** :math:`L_a` **from equation** :eq:`LA` **)**

.. math:: {{(\Delta y_j)^2} \over {2 tan \alpha}} \le  f_a V_a + f_b V_b
          \le  {{(\Delta y_j + \Delta y_{j+1})^2} \over {2 tan \alpha}}

Area of a triangle:

.. math:: L_w = \left ({2 f_b V_b} \over {tan \alpha}  \right )^{0.5}
   :label: LWB

.. image:: ../_images/limits_b.png
   :scale: 75%
   :align: center

**Case C solution for** :math:`L_w` **(** :math:`L_a` **from equation** :eq:`LA` **)**

.. math::  {{(\Delta y_j + \Delta y_{j+1})^2} \over {2 tan \alpha}} \le  f_a V_a + f_b V_b
          \le (\Delta y_{j+1} + \Delta y_{j}) \Delta x- {{(\Delta y_j + \Delta y_{j+1})^2} \over {2 tan \alpha}}

From the area of a trapezoid:

.. math:: L_w = {{f_a V_a + f_b V_b} \over {\Delta y_j + \Delta y_{j+1}}}
   :label: LWC

.. image:: ../_images/limits_c.png
   :scale: 75%
   :align: center

**Case D solution for** :math:`L_w` **(this is just Case B but rotated 180 degrees and i.t.o. air** :math:`L_a` **from equation** :eq:`LA` **)**

From Case B (i.t.o. air):

.. math:: {{(\Delta y_{j+1})^2} \over {2 tan \alpha}} \le  (1-f_a) V_a + (1-f_b) V_b
          \le  {{(\Delta y_j + \Delta y_{j+1})^2} \over {2 tan \alpha}}

From the area of a triangle:

.. math:: L_w = \Delta x - \left ({2 (1-f_a) V_a} \over {tan \alpha}  \right )^{0.5}
   :label: LWD

.. image:: ../_images/limits_d.png
   :scale: 75%
   :align: center

**Case E  solution for** :math:`L_w` **(this is just Case A but rotated 180 degrees and i.t.o. air** :math:`L_a` **from equation** :eq:`LA` **)**

From Case A (i.t.o. air):

.. math:: (1-f_a) V_a + (1-f_b) V_b \le  {({\Delta y_{j+1})^2} \over {2 tan \alpha}}

From the area of a triangle:

.. math:: L_w = \Delta x
   :label: LWE

.. math:: f_a = 1

.. image:: ../_images/limits_e.png
   :scale: 75%
   :align: center

:math:`L_w` and :math:`L_a` on south, east and west faces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Solution for** :math:`L_w` **on south face by replacing** :math:`j` **by** :math:`j-1` **in equations** :eq:`FA`, :eq:`FB`, :eq:`LW`, :eq:`LA`, :eq:`LWA`, :eq:`LWB`, :eq:`LWC`, :eq:`LWD`, :eq:`LWE`

* **Solution for** :math:`L_w` **on east face**:

.. math::  L_w = f_{i+1,j} \Delta y_j
   :label: LWEast

* **Solution for** :math:`L_w` **on west face by replacing** :math:`i` **by** :math:`i-1` **in** :eq:`LWEast`

* **Solution for** :math:`L_a` **on south face from equation** :eq:`LA`

* **Solution for** :math:`L_a` **on west and east face**

.. math:: L_a = \Delta y_j - L_w

:math:`\mu` on north face (:math:`\mu_n`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**What do we need to compute?**

The diffusive flux on the north face:

.. math::
  
   D_n = \mu_n {{A_n} \over {\Delta y_n}}

What is :math:`\mu_n`?

* We are concerned with the viscosity at the interface between points P and N.
* As Patankar (1980) suggests, we must determine the distance to the interface from the points P and N.
* **The governing assumption here is that very steep angles > 45 degrees or < -45 degrees are not present**
* Note that this assumption only applies to the viscosity and not to the convective flux, which does allow sharp angles
* **The second governing assumption is that the viscosity at the interface is represented by the harmonic mean, so significant variation in the viscosity (e.g. due to turbulence) between the grid points at the interface is negligible**

The use of the harmonic mean follows the suggestion of Patankar (1980) who used it for non-uniform conductivity in the diffusive term. Patankar (1980) dismisses the use of the arithmetic mean as too simplistic, but he did not consider the geometric mean or the infinite norm mean see Schmeling et al (2008). We also don't know whether turbulent viscosity can be best represented by these averaging methods. 

There are four cases:

* :math:`f_a < 0.5` and :math:`f_b = 0` so :math:`d_w = 0`
* :math:`f_a > 0.5` and :math:`f_b = 0` so :math:`d_w = (f_a - 0.5)\Delta y_j`
* :math:`f_a > 0.5` and :math:`f_b < 0.5` so :math:`d_w = {(f_a - 0.5)\Delta y_j} +f_b \Delta y_{j+1}`
* :math:`f_a > 0.5` and :math:`f_b > 0.5` so :math:`d_w = {(f_a - 0.5)\Delta y_j} +0.5 \Delta y_{j+1}`
 
This can be summarized like this, where :math:`d_w` is the distance from P to the interface:

.. math:: d_w = max(f_a-0.5,0)\Delta y_j + min(f_b,0.5)\Delta y_{j+1}
   :label: DW

.. image:: ../_images/viscosity.png
   :scale: 75%
   :align: center

The distance from N to the interface is :math:`d_a` and this is simply the distance between P and N minus :math:`d_w`:

.. math:: d_a = {{\Delta y_j + \Delta y_{j+1}} \over 2} - d_w
   :label: DA

**Harmonic Mean Viscosity**

This is how Patankar (1980) defines the harmonic mean viscosity (for the north face):

.. math:: \mu_n = \left ( {f_w \over \mu_w} + {f_a \over \mu_a} \right )^{-1}
   :label: MUN

Where:

.. math:: f_w = {d_w \over {d_a + d_w}}
   :label: FWMU

And:

.. math:: f_a = {d_a \over {d_a + d_w}}
   :label: FAMU


Substituting Equation :eq:`FWMU` and :eq:`FAMU` into :eq:`MUN` and re-arranging gives:

.. math:: \mu_n ={ {\mu_w \mu_a (d_w + d_a)} \over {\mu_a d_w + \mu_w d_a}}
   :label: MUNORTH

:math:`\mu` on the south, east and west faces (:math:`\mu_s`, :math:`\mu_e`, :math:`\mu_w`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Solution for** :math:`\mu_s` **on south face by replacing** :math:`j` **by** :math:`j-1` **in equation** :eq:`DW` **and** :eq:`DA`

* **Solution for** :math:`\mu_e` **on east face by convectional VOF method (the arithmetic mean) assuming linear variations in viscosity**:

.. math:: \mu_e = f_{i+1,j} \mu_w + (1-f_{i+1,j}) \mu_a
   :label: MUE

* **Solution for** :math:`\mu_w` **on west face by replacing** :math:`i` **by** :math:`i-1` **in equation** :eq:`MUE`


Total mass in u-control volume
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total mass in the u-control volume can be seen as a function of time, because the values of :math:`f` are changing with time, especially near the interface. So the mass in the control volume is:

.. math:: m = \int_{V} \rho dV = (f_a \rho_w + f_b \rho_a)V

.. math:: f_b = 1-f_a

(Seriously consider changing :math:`f_a` to :math:`f_P` and :math:`f_b` to :math:`f_N` to avoid confusion with air)

v-control volume
----------------

The angle :math:`\alpha` is from Brackbill et al (1994) - **I'm not sure why the angle the interface makes with the x-axis would be different for the v control volume? Or why we would need another equation for the same value? Is this because the v-volume fraction is stored at the faces and not the cell centres?**

.. math::

  tan(\alpha) = \left \vert {{(f_{i+1,j+2}+2f_{i+1,j+1}+f_{i+1,j}-f_{i-1,j+2}-2f_{i-1,j+1}-f_{i-1,j})\Delta y} \over
  {(f_{i+1,j+2}+2f_{i,j+2}+f_{i-1,j+2}-f_{i+1,j}-2f_{i,j}-f_{i-1,j})\Delta x}} \right \vert

**Why isn't this the same as Equation** :eq:`alpha`?

:math:`L_w` and :math:`L_a` where :math:`tan(\alpha) < \Delta y / \Delta x` on north face
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Case A - Solution for** :math:`L_w`

.. math:: f_{i,j+1} \le {{\Delta x tan \alpha} \over {2 \Delta y}}

.. math:: L_w = 0


**Case B - Solution for** :math:`L_w`

.. math::  {{\Delta x tan \alpha} \over {2 \Delta y}} \le f_{i,j+1} \le  {1 \over 2} +{{\Delta x tan \alpha} \over {2 \Delta y}}

From the area of a trapezoid:

.. math:: L_w = \left ( {{f_{i,j+1}\Delta y } \over {\Delta x tan \alpha} } - {1 \over 2}   \right )


**Case C - Solution for** :math:`L_w`

.. math:: f_{i,j+1} \ge {1 \over 2} +{{\Delta x tan \alpha} \over {2 \Delta y}}

.. math:: L_w = \Delta x

.. image:: ../_images/v_con_1.png
   :scale: 75%
   :align: center

:math:`L_w` and :math:`L_a` where :math:`tan(\alpha) > \Delta y / \Delta x` on north face
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Case A - Solution for** :math:`L_w`

.. math:: f_{i,j+1} \le {{\Delta y} \over {8 \Delta x tan \alpha}}

.. math:: L_w = 0


**Case B - Solution for** :math:`L_w`

.. math::  {{\Delta y} \over {8 \Delta x tan \alpha}} \le f_{i,j+1} \le {{\Delta y} \over {2 \Delta x tan \alpha}}

From the area of a triangle:

.. math:: L_w ={{(2 \Delta x \Delta y f_{i,j+1} tan \alpha)^{0.5} - {\Delta y \over 2}} \over {tan \alpha}}


**Case C - Solution for** :math:`L_w`

.. math:: {{\Delta y} \over {2 \Delta x tan \alpha}} \le f_{i,j+1} \le 1-{{\Delta y} \over {2 \Delta x tan \alpha}} 

.. math:: L_w =  f_{i,j+1} \Delta x

**Case D - Solution for** :math:`L_w`

.. math:: 1-{{\Delta y} \over {2 \Delta x tan \alpha}} \le f_{i,j+1} \le 1-{{\Delta y} \over {8 \Delta x tan \alpha}} 

.. math:: L_w = \Delta x + {\Delta y \over {2 tan \alpha}} - {{(2 \Delta x \Delta y (1-f_{i,j+1}) tan \alpha)^{0.5}} \over {tan \alpha}}

**Case E - Solution for** :math:`L_w`

.. math:: f_{i,j+1} \ge 1-{{\Delta y} \over {8 \Delta x tan \alpha}}

.. math:: L_w = \Delta x

.. image:: ../_images/v_con_2.png
   :scale: 75%
   :align: center


**Resulting flux**

Depending on the value of :math:`f_{i,j+1}` and :math:`\alpha`, the fluid flux produced may vary from:

* :math:`\rho_a v_n \Delta x` when :math:`L_w = 0` and :math:`L_a = \Delta x` 
* :math:`\rho_w v_n \Delta x` when :math:`L_w = \Delta x` and :math:`L_a = 0` 


:math:`L_w` and :math:`L_a` on south, east and west faces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Solution for** :math:`L_w` **on south face by replacing** :math:`j` **by** :math:`j-1` **in equations for north face**

* **Solution for** :math:`L_w` **on east face**

.. math:: f_a = {f_{i,j}+f_{i+1,j}}\over 2

.. math:: f_b = {f_{i,j+1}+f_{i+1,j+1}}\over 2

.. image:: ../_images/v_con_east.png
   :scale: 75%
   :align: center

.. math:: L_w = max(f_a-0.5,0)\Delta y_j+min(f_b,0.5)\Delta y_{j+1}

.. image:: ../_images/v_con_west.png
   :scale: 75%
   :align: center

**not sure about the following:**

.. math:: L_e ={{{ {{\Delta y_{j+1}}}\over 2} + {{\Delta y_{j}}\over 2}} - L_w}

* **Solution for** :math:`L_w` **on west face by replacing** :math:`i-1` **with** :math:`i` **in equations for east face**


:math:`\mu` on north face (:math:`\mu_n`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: d_w = f_{i,j+1} {\Delta y_{j+1}}
   :label: DWN

.. math:: d_a = {\Delta y_{j+1}} - d_w
   :label: DAN

**I'm not sure about the positioning of the indices, in the above equations**

**BIG QUESTION: Why is the wet length based on the u-volume fraction cell-centred and the wet length based on the v-volume fraction face-centred?**

.. image:: ../_images/v_con_4.png
   :scale: 75%
   :align: center

The harmonic mean is used for the viscosity :math:`\mu_n` from Equation :eq:`MUNORTH`

:math:`\mu` on the south, east and west faces (:math:`\mu_s`, :math:`\mu_e`, :math:`\mu_w`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Solution for** :math:`\mu_s` **on south face by replacing** :math:`j+1` **by** :math:`j` **in equation** :eq:`DWN` **and** :eq:`DAN`

* **Solution for** :math:`\mu_e` **on east face by convectional VOF method (the arithmetic mean) assuming linear variations in viscosity**:

.. math:: \mu_e = f_e \mu_w + (1-f_e) \mu_a

Similar to :math:`L_w`:

.. math:: f_e = max(f_a-0.5) + min(f_b, 0.5)


* **Solution for** :math:`\mu_w` **on west face by replacing** :math:`i` **by** :math:`i-1` **in equations for east face**


Total mass in v-control volume
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The total mass in the v-control volume:

.. math:: m = \int_{V} \rho dV = (f \rho_w + (1-f) \rho_a)V

Similar to :math:`L_w`:

.. math:: f = max(f_{i,j}-0.5) + min(f_{i,j+1}, 0.5)

It should be mentioned that there are various methods for reconstructing the interface, e.g. Rider and Kothe (1998)

Solution Algorithms
===================

SIMPLE and SIMPLEC Algorithms
-----------------------------

* The u momentum equation :eq:`u_momentum_equation` and v momentum equation :eq:`v_momentum_equation` are solved iteratively until the normalised residuals are small (perhaps using the tri-diagonal matrix algorithm). At this point we have satisfied momentum conservation.
* However, there is no guarantee that this solution satisfies the **volume conservation equation**, from equation :eq:`volume`:

.. math:: u_e A_e - u_w A_w + v_n A_n v_s A_s = 0

* A pressure-correction algorithm is needed in the loop in order to satisfy volume conservation by correcting the following according to the volume conservation equation:

  - velocity
  - pressure

SIMPLE
~~~~~~

The SIMPLE algorithm leads to a slow convergence when there is a sudden change in the geometry or solution property. The u velocity correction is :math:`u_P^c`:

.. math:: u_P^c ={ {p_w^c A_w - p_e^c A_e} \over {{m_P^n \over \Delta t} + \sum a_N}}

SIMPLEC
~~~~~~~

The SIMPLEC algorithm is used instead:

.. math:: u_P^c ={ {p_w^c A_w - p_e^c A_e} \over {{m_P^n \over \Delta t}}}

Now obtain the pressure correction by replacing :math:`u = u^* + u^c` (:math:`u` is the true velocity and equals the uncorrected velocity :math:`u^*` plus the corrected velocity :math:`u^c`)

.. math:: p_P^c = {1 \over a_P}{(a_E p_E^c + a_W p_W^c + a_N p_N^c + a_S p_S^c + b)}

where b is the mass imbalance around P:

.. math:: b = -(u_e^*A_e - u_w^*A_w + v_n^*A_n - v_s^*A_s)

The pressure and velocity corrections are applied like this:

.. math:: p_P = p_P^* + p_P^c

.. math:: u_P = u_P^* + u_P^c

.. math:: v_P = v_P^* + v_P^c

 
CICSAM Algorithm
----------------

* At this point we have satisfied the equations for volume, mass, u momentum and v momentum conservation.

* The final equation to satisfy is the equation for volume fraction conservation :eq:`volume_fraction`, which is solved explicitly:

.. math:: f_P = f_P^n - {\Delta t \over V}(u_e^n f_e^n A_e - u_w^n f_w^n A_w + u_n^n f_n^n A_n - u_s^n f_s^n A_s)
   :label: volume_fraction_solution

* The values on the surfaces are computed by a high resolution scheme CICSAM. The reason for using CICSAM is that this method is able to produce a sharp water surface

* Normal to the interface, the volume fraction is between 0 and 1 on only one node.

* The problem with Equation :eq:`volume_fraction_solution` is that when the points :math:`u, w, n, s` and :math:`P` are all in water, the RHS of Equation :eq:`volume_fraction_solution` is:

.. math:: f_P = 1 - {\Delta t \over V}(u_e^n A_e - u_w^n A_w + u_n^n A_n - u_s^n A_s)

* Clearly :math:`f_P` won't be exactly 1, so the equation for volume fraction conservation won't reduce to the equation of volume conservation (which it should in this case). So in solving Equation :eq:`volume_fraction_solution` a correction is applied to the velocity:

.. math:: V^c = {{(u_e^n A_e - u_w^n A_w + v_n^n A_n - v_s^n A_s)} \over {2 \Delta x + 2 \Delta y}} 

* The correction is applied like this to Equation :eq:`volume_fraction_solution`

.. math:: u_e^n = u_e^n + V^c
.. math:: u_w^n = u_w^n - V^c
.. math:: u_n^n = u_n^n + V^c
.. math:: u_s^n = u_s^n - V^c

CICSAM, SIMPLE-PISO Procedure
-----------------------------

In order to accelerate convergence, we add a second velocity and pressure produced by the PISO algorithm **at the last iteration in each time step**

The implicit SIMPLEC-PISO algorithm is combined with the explicit CICSAM algorithm like this:

1) The solution of volume fraction is calculated from the solutions of :math:`u^n`, :math:`v^n` and :math:`f^n` by CICSAM.
2) The solutions of :math:`u`,  :math:`v` and :math:`p` are obtained from the SIMPLEC-PISO algorithm after :math:`f` has already been obtained.

Application to a Non-Linear, Non-Breaking Deep Water Wave
=========================================================

* A wave is modelled using the inviscid non-linear solution to the wave equation as the inlet and initial conditions (the harmonic solution to the velocity potential described by Laplace's equation plus boundary conditions).
* The model uses :math:`h/L = 0.6` i.e. it's a deep water case (:math:`h/L > 0.5`)
* Also :math:`ka = 0.1885` i.e. it's non-linear (:math:`ka > 0.1`)
* And :math:`2a/L = 0.06` i.e. it's non-breaking (:math:`2a/L < 0.08`)
* No wind input is used, i.e. the air follows the water

* The model has a wall at the bottom boundary, a symmetry condition at the top, an inlet velocity, plus an opening at the outlet.

Velocity
--------

* Potential flow compares well with the wet-dry areas method in the water, however the difference occurs in the air and near the interface:
   - The velocity field is discontinuous in the potential flow solution, because the viscosity is zero
   - The wet-dry areas method produces a rapid but continuous change in the velocity field
* However remarkably over most of the velocity profile the difference is less than 1% 
* The v-velocity variation is continuous in both cases and little difference is seen

Streamlines
-----------

* At the peak and trough, streamlines for potential flow are centred at the interface, whereas streamlines for the wet-dry areas method are centred in the air

**Generally speaking this exercise of comparing the potential flow solution with the wet-dry areas method was for validation purposes and it showed that sensible results are being produced**

**I wonder whether a comparison between the standard arithmetic approach to the mass fluxes and the wet-dry areas method would also be useful to compare the shear stress computed at the interface.**

Application to Collapse of a Water Column
=========================================

* The validation with experimental data was also conducted against images for the a collapse of a water column. 
* It was found that the qualitative flow features, such as the hump in the water column, the strong jetting and the air pocket near the wall were all present in both the simulation and the experimental images. However, as the paper admits, there is no quantitative measurement to make this validation more complete.

Application to a Non-Linear, Breaking Shallow Water Wave
========================================================

* The breaking point is defined as the point where the speed of the water particles becomes larger than the phase speed of the wave.
* At this point the front of the water wave becomes vertical.
* To test whether the physics of this phenomenon is reproduced, a numerical simulation was run with the wet-dry areas method using the same parameters as the PIV experiment of Change and Liu (1998), which is where wave breaking was measured.
* Wave steepness: :math:`2a/L = 0.12` (i.e. :math:`2a/L > 0.08` which implies a breaking wave)
* Depth: :math:`h/L = 0.166` (i.e. :math:`h/L < 0.5` which implies a shallow wave)
* It was found that the dimensionless water particle velocity was in the range 1.0c to 1.1c (where c is the phase speed), which matches the measured value of 1.07c.

Grid Independence, Robustness and Efficiency of the Wet-Dry Areas Method
========================================================================

Grid Independence
-----------------

* This was performed for the water column collapse application for T=3.99. 
* Uniform, structured grids were used in this case.
* **I wonder whether non-uniform structured grids or non-uniform unstructured grids might be more numerically efficient, although they may not be less diffuse.**

The Robustness
--------------

* The first 20 timesteps of the water column collapse application was also used test robustness.
* The SIMPLEC algorithm used is unconditionally stable, however the convergence at each timestep depends on the timestep - i.e. divergence may occur if the timestep is too large.
* The robustness was tested by using a fixed timestep and testing to see if the volume residual became less than :math:`5x10^{-4}`
* The minimum number of iterations per timestep was set to 15, which was achieved quickly after the first five timesteps.

Efficiency
----------

* The efficiency of the wet-dry areas method is less than the average density method due to the use of different cases in the u and v control volumes
* For the water column collapse it is 20% greater (as this uses all 4 cases)
* For the non-linear non-breaking wave and the non-linear breaking wave it is around 10% longer (as these are mainly case 1)

Conclusions
===========

The Navier-Stokes Equations have been solved using:

* Conservative integral form
* Standard staggered grid
* A new analytical expression for wet-dry areas in the convective mass flux (instead of the standard average density approach), which is a function of volume fraction
* SIMPLEC algorithm (which is faster than SIMPLE for interfacial flows)

Advantages of wet-dry areas method:

* Improves the computation of convective mass flux
* Bussmann's consistency problem is avoided
* The new numerical method is robust 

Disadvantages of wet-dry areas method:

* Less efficient than the average density approach.

