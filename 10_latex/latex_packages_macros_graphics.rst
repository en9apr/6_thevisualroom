=============================
Packages, Macros and Graphics
=============================

.. contents::
   :local:

.. highlight:: latex

Packages
========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			$\usepackage{fullpage}$
			 
     - Full page (minimal margins) - positioned after ``\documentclass`` and before ``\begin{document}`` 
   * - ::

			$\usepackage[top=1in, bottom=1in, left=1in, right=1in]{geometry}$

     - Specification of the margins using the geometry package, same as ``\usepackage{fullpage}``.
   * - ::

                        $\usepackage[margin=1in, paperwidth=8.5in, paperheight=11in]{geometry}$

     - Specification of margins using the geometry package with specification of paper width and height, same as ``\usepackage{fullpage}``.
	   
   * - ::

			$\usepackage{amsfonts}$

     - Used for displaying symbols such as the natural number symbol (``$\mathbb{N}$``), the real number symbol(``$\mathbb{R}$``) and the integer symbol (``$\mathbb{Z}$``).


Macros
======

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			\def\eq1{y=\frac{x}{3x^2+x+1}}
     - ``\eq1`` is the name of the equation (can be anything, but not a preset name). The use of this equation is ``$\eq1$``. The macros are positioned in the preamble after ``\documentclass`` and before ``\begin{document}``.  
   * - ::

			\def\labelaxes{Remember, remember the fifth of November}
     - ``\labelaxes`` is the name of the string (can be anything, but not a preset name). The use of this string is ``$\labelaxes$``. The macros are positioned in the preamble after ``\documentclass`` and before ``\begin{document}``. 

Graphics
========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			\usepackage{graphicx}

     -  Can also centre the image  ``\includegraphics[width=5in]{quad1.png}``

   * - ::

			\includegraphics{quad1.png}

     - Use of graphix package to display an image. Will look for the file in the same directory as the .tex file. Can only use .png, .jpg, .gif or .pdf files (not sure - we can also use .eps?)

   * - ::

			\includegraphics[width=5in]{quad1.png}

     - Use of graphix package to display an image, specifying the width.

   * - ::

                        \begin{center}
			    \includegraphics{quad1.png}
			\end{center}

     - Use of graphix package to display an image, centered.

   * - ::
                       
		        \includegraphics[scale=0.5]{quad1.png}
		

     - Use of graphix package to display an image, scaled to 50%.

   * - ::
                       
		        \includegraphics[angle=45]{quad1.png}
		

     - Use of graphix package to display an image, rotated by 45 degrees.
