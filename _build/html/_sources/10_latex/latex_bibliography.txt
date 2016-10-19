==================
 LaTeX Bibiography
==================

This is an example of a LaTeX bibliography.

.. contents::
   :local:

.. highlight:: latex

::

   %++++++++++++++++++++++++++++++++++++++++
   % References section will be created automatically 
   % with inclusion of "thebibliography" environment
   % as it shown below. See text starting with line
   % \begin{thebibliography}{99}
   % Note: with this approach it is YOUR responsibility to put them in order
   % of appearance.

   \begin{thebibliography}{99}

   \bibitem{melissinos}
   A.~C. Melissinos and J. Napolitano, \textit{Experiments in Modern Physics},
   (Academic Press, New York, 2003).

   \bibitem{Cyr}
   N.\ Cyr, M.\ T$\hat{e}$tu, and M.\ Breton,
   % "All-optical microwave frequency standard: a proposal,"
   IEEE Trans.\ Instrum.\ Meas.\ \textbf{42}, 640 (1993).

   \bibitem{Wiki} \emph{Expected value},  available at
   \texttt{http://en.wikipedia.org/wiki/Expected\_value}.

   \end{thebibliography}

The references are cited like this (~ is a space):

::

   ~\cite{melissinos, Cyr, Wiki}


