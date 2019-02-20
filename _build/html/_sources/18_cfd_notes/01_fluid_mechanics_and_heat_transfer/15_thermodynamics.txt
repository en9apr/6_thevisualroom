==============
Thermodynamics
==============

.. contents::
   :local:

.. highlight:: latex

What are examples of path functions and state variables?
========================================================

* Path functions: Heat and work (about the boundary)
* State variables: Internal energy, enthalpy, entropy (about the fluid)

What is the First Law of Thermodynamics and what kind of processes is it applicable to?
=======================================================================================

.. math::

    \text{Net heat added to system} + \text{Net work done on system} = \text{Increase in internal energy of system}
    
.. math::

    dq + dw = de
    
(lower case implies per unit mass)

The first law applies to reversible and irreversible processes

What is the Second Law of Thermodynamics and what is it for a reversible process?
=================================================================================

.. math::

    \text{Gross heat added to system} + \text{Net work done by system}  \ge 0
    
.. math::

    dq + \Sigma = Tds

i.e. there is irreversible work done by internal molecular changes :math:`\Sigma`

For a reversible process :math:`\Sigma = 0`

i.e. :math:`ds = {dq \over T}`

What is Enthalpy?
=================

From the :math:`1^{st}` law, for both reversible and irreversible processes:

.. math::

    dq = d(e+pv) = dh
    
where :math:`dw = -d(pv)`

Enthalpy is a state variables

For both reversible and irreversible processes it equals the internal energy plus pressure-volume potential energy.

What is Entropy?
================

From the :math:`2^{nd}` law for a reversible process:

.. math::

    s_2 - s_1 = \int_1^2 {dq \over T}
    
Entropy is a macro state variable

For a reversible process it is the infinitisimal heat added over constant temperature (non-adiabatic).

What are the two property relations for entropy change from the 1st and 2nd laws of thermodynamics at constant pressure and constant volume?
============================================================================================================================================

* 1st Law (reversible and irreversible) :math:`\longrightarrow de = dw+dq`
* Reversible work done by system :math:`\longrightarrow dw = -pdv`
* Reversible heat added to system (2nd law) :math:`\longrightarrow dq = Tds`
* For reversible and irreversible processes (1) :math:`\longrightarrow Tds = de + pdv`
* Definition of enthalpy :math:`\longrightarrow dh = d(e+pv) = de + d(pv) = de + pdv + vdp`
* By rearrangement :math:`\longrightarrow de = dh - pdv - vdp`
* For reversible and irreversible processes (2) :math:`\longrightarrow Tds = dh - vdp`

What is the entropy change for an ideal gas for a reversible and irreversible process?
======================================================================================

* :math:`ds = {de \over T}+{p \over T}dv`
* :math:`ds = {dh \over T} - {v \over T}dp`

* Calorically perfect 

:math:`\longrightarrow de=c_vdT`

:math:`\longrightarrow dh=c_pdT`

* Ideal gas

:math:`\longrightarrow p = {RT \over v}`

:math:`\longrightarrow v = {RT \over p}`

* Hence

:math:`\longrightarrow s_2 - s_1 = c_v \int_1^2 {dT \over T} + R \int_1^2 {dv \over v} = c_v ln{T_2 \over T_1} + R ln{v_2 \over v_1} = c_v ln{T_2 \over T_1} + R ln{\rho_1 \over \rho_2}`

:math:`\longrightarrow s_2 - s_1 = c_p \int_1^2 {dT \over T} - R \int_1^2 {dp \over p} = c_p ln{T_2 \over T_1} - R ln{p_2 \over p_1}`

How are the isentropic relations for an ideal gas derived from the entropy change?
==================================================================================

* Isentropic :math:`\longrightarrow s_2 - s_1 = 0`

* Hence

:math:`\longrightarrow 0 = c_v ln{T_2 \over T_1} - R ln{\rho_2 \over \rho_1}`

:math:`\longrightarrow 0 = c_p ln{T_2 \over T_1} - R ln{p_2 \over p_1}`

* Definition of :math:`c_v` and :math:`c_p`:

:math:`c_v = {R \over {\gamma -1}}`

:math:`c_v = {{\gamma R} \over {\gamma -1}}`

* Hence:

:math:`\longrightarrow ln{\rho_2 \over \rho_1} = {1 \over {\gamma -1}} ln{T_2 \over T_1}` 

:math:`\longrightarrow ln{p_2 \over p_1} = {{\gamma} \over {\gamma -1}} ln{T_2 \over T_1}`

* Hence:

:math:`\longrightarrow {\rho_2 \over \rho_1} =  {T_2 \over T_1}^{1 \over {\gamma -1}}` 

:math:`\longrightarrow {p_2 \over p_1} =  {T_2 \over T_1}^{{\gamma} \over {\gamma -1}}`

* Hence:

:math:`\longrightarrow {p_2 \over p_1} =  {\rho_1 \over \rho_2}^{\gamma}`

:math:`\longrightarrow {p \over {\rho^{\gamma}}} =  \alpha` 







