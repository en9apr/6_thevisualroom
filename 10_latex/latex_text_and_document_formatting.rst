============================
Text and Document Formatting
============================

.. contents::
   :local:

.. highlight:: latex

Text format
===========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

        This will produce \textit{italicized} text.

        This will produce \textbf{bold-faced} text.

        This will produce \textsc{small caps} text.

        Please visit Mrs. Krummel's website at \texttt{http://mrskrummel.com}.
			 
     - 
       * Italicized text
       * Bold faced text
       * Small caps text
       * Typewriter font

Text size
=========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

        Please excuse my \begin{tiny}dear Aunt Sally\end{tiny}.

        Please excuse my \begin{small}dear Aunt Sally\end{small}.

        Please excuse my dear Aunt Sally.

        Please excuse my \begin{large}dear Aunt Sally\end{large}.

        Please excuse my \begin{Large}dear Aunt Sally\end{Large}.

        Please excuse my \begin{huge}dear Aunt Sally\end{huge}.

        Please excuse my \begin{Huge}dear Aunt Sally\end{Huge}.

     -        
       * tiny
       * small
       * normal
       * large
       * Large
       * huge
       * Huge

Text alignment
==============

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

        \begin{center} This is centered \end{center}

        \begin{flushleft} This is left-justified. \end{flushleft}

        \begin{flushright} This is right-justified. \end{flushright}

     - 
       * Centred
       * Left
       * Right

Title
=====

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::
       
        \title{My Practice \LaTeX \ Document}
	\author{Andrew Roberts}
	\date{\today}
	\maketitle

     - 
       * Title:

         - LaTeX has a special display style ``\LaTeX``
	 - ``\`` is a space

       * Author
       * Date (optional)
       * ``\maketitle`` actually creates the title


Sections
========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::
       
        \section{Linear Functions}
	    \subsection{Slope-Intercept Form}
	    The slope-intercept form of a linear function is given by $y=ax+b$.
            \subsection{Standard Form}
            \subsection{Point-Slope Form}
	\section{Quadratic Functions}
            \subsection{Vertex Form}
            \subsection{Standard Form}
            \subsection{Factored Form}
     - 
       * Section 1
       * Subsection 1.1
       * Subsection 1.2
       * Subsection 1.3
       * Section 2
       * Subsection 2.1
       * Subsection 2.2
       * Subsection 2.3

Table of Contents
=================

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::
       
        \tableofcontents

     - Creates a table of contents. Note that you have to quick build twice in Texmaker to get this to work. One click will collect the information, the second click will actually create the table. 
