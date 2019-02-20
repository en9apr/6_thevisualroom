=====
Lists
=====

.. contents::
   :local:

.. highlight:: latex

Enumerated
==========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			\begin{enumerate}
			   \item pencil
			   \item paper
			   \item calculator
			   \item ruler
			   \item notebook
			      \begin{enumerate}
			         \item assessments
			            \begin{enumerate}
				       \item tests
				       \item quizzes
				    \end{enumerate}
			         \item homework
			         \item notes
			      \end{enumerate}
			   \item graph paper
			\end{enumerate}
			 
     - Enumerated list:

       * Numbered list (first level)
       * Alphabetical list (second level)
       * Roman numeral list (third level)
	   

Itemized
========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

			\begin{itemize}
			   \item pencil
			   \item paper
			   \item calculator
			   \item ruler
			   \item notebook
			      \begin{itemize}
			         \item assessments
			            \begin{itemize}
				       \item tests
				       \item quizzes
				    \end{itemize}
			         \item homework
			         \item notes
			      \end{itemize}
			   \item graph paper
			\end{itemize}

     - Itemized list:
       
       * Round bullets (first level)
       * Minus sign (second level)
       * Star (third level)

Labelled
========

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

	\begin{enumerate}
           \item[Commutative] $a+b=b+a$
           \item[Associative] $a+(b+c)=(a+b)+c$
           \item[Distributive] $a(b+c)=ab + ac$
        \end{enumerate}

     - Labelled list (name in square brackets) - right justified labels
