=========================
Creating a LaTeX Document
=========================

.. contents::
   :local:

.. highlight:: latex

Basic Document Setup
====================

The document is setup like this:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			\documentclass[11pt]{article}
		   
     - The preamble
   * - ::

			\begin{document}
			This is the body of the document
			\end{document}

     - Body of the document
   * - ::

			Inline maths $(x+1)$

     - Inline maths
   * - ::

			Display maths $$(x+1)$$

     - Display (centred on it's own line) maths
   * - ::

			Line 1. \\
			Line 2.

     - New line

Files
=====

There are five files created, the two most important ones are in bold (.pdf and .tex)

* filename.aux
* filename.txt
* **filename.pdf** (output) 
* filename.synctex.gz
* **filename.tex** (latex command file)

Templates
=========

Useful report template for LaTeX:

`Overleaf Sample Lab Report for PHYS 349 at the University of Redlands <https://www.overleaf.com/latex/templates/sample-lab-report-for-u-of-r-phys-349/pgsyqngcyjxk#.VhQxeXpVhBd>`_







