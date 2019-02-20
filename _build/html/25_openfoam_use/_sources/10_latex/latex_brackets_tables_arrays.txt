===========================
Brackets, Tables and Arrays
===========================

.. contents::
   :local:

.. highlight:: latex

Brackets
========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			$(x+1)$
			 
     - Brackets
	   
   * - ::

			$3[2+(x+1)]$

     - Square brackets and normal brackets
   * - ::

			$\{a,b,c\}$

     - Curly brackets - **curly bracket is reserved in Latex**
   * - ::

			$\% 50$

     - Dollar sign - **dollar sign is reserved in Latex**
   * - ::

			$3\left(\frac{2}{5}\right)$

     - Expanded brackets
   * - ::

			$3\left[\frac{2}{5}\right]$

     - Expanded square brackets
   * - ::

			$3\left\{\frac{2}{5}\right\}$

     - Expanded curly brackets
   * - ::

			$\left|\frac{2}{5}\right|$

     - Absolute value
   * - ::

			$\left.\frac{2}{5}\right|$

     - No left bracket
   * - ::

			$\left.\frac{dy}{dx}\right|_{x=1}$

     - Value of a derivative at x=1

Tables
======

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			\begin{tabular}{|c|c|c|c|c|c|}
			\hline
			$x$ & 1 & 2 & 3 & 4 & 5 \\ \hline
			$f(x)$ & 10 & 11 & 12 & 13 & 14 \\ \hline
			\end{tabular}
     - Table with border and grid:

       * ``c`` = column
       * ``&`` = next column
       * ``|`` = vertical line
       * ``\\`` = new line
       * ``\hline`` = horizontal line

Equation Arrays
===============

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			\begin{eqnarray*}
			5x^2 - 9 &=& x + 3 \\
			4x^2 &=& 12 \\
			x^2 &=& 3 \\
			x &\approx& \pm 1.732
			\end{eqnarray*}

     - Equation array:

       * ``*`` = no equation numbers
       * ``\\`` = new line
       * ``&=&`` = line up by equals sign
       * ``&\approx&`` = line up by approx. equals sign
       * ``\pm`` = plus or minus

