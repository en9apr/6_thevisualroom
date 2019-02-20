===============
Continuum Model
===============

.. contents::
   :local:

.. highlight:: latex

What is the Continuum Model?
============================

Discontinuous molecular collisions are replaced with continuous statistically averaged values, e.g.

* Pressure = Average collision force per unit area.

Why is the Continuum Model used?
================================

Average fluid properties in the Continuum model are used because:

* Molecular density of fluid is large with respect to characteristic length.
* Interaction forces (potentials) are now known.

What does the Continuum Model enable?
=====================================

Enables the conservation laws:

* Mass
* Momentum
* Energy

Where can the Continuum Model be applied?
=========================================

Can be applied to:

* Finite volume element
* Infinitesimal differential element

What is the Knudsen number?
===========================

.. math:: {Kn = {\lambda \over L}}

where:

* :math:`\lambda` = mean free path
* :math:`L` = characteristic length (e.g. pipe diameter, boundary layer thickness, shockwave thickness)

What is the Mean Free Path?
===========================

.. math:: {\lambda = {1 \over {\sqrt{2} n S}}}

where:

* :math:`n` = number of particles per :math:`m^3` from the Boltzmann equation (:math:`p=nk_bT`)
* :math:`S` = surface area of a spherical molecule = :math:`\pi d^2`
* :math:`d` = molecular diameter or collision diameter

What are the limitations of the Continuum Model?
================================================

The Continuum Model is valid for Knudsen Numbers :math:`\ll 1` i.e.

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Knudsen Number
     - Meaning
   * - :math:`Kn > 1`
     - Molecular dynamics (e.g. gas in nano flow)
   * - :math:`0.1 < Kn < 1`
     - Transition/continuum mechanics
   * - :math:`Kn \sim 0.01`
     - Continuum (NS) + slip boundary
   * - :math:`Kn \ll 1`
     - Continuum (NS)
   * - :math:`Kn \to 0`
     - Continuum (Euler)   
     
What does the Continuum Model not permit?
=========================================

Mass-energy conversion e.g.

* Nuclear reactions
* Relativity