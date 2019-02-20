================
 LaTeX Preamble
================

This is an example of a LaTeX pre-amble set for A4 paper with 1 inch margins and 12pt font.

.. contents::
   :local:

.. highlight:: latex

::


   %++++++++++++++++++++++++++++++++++++++++
   % Don't modify this section unless you know what you're doing!
   \documentclass[a4paper,12pt]{article}
   \usepackage{tabularx} % extra features for tabular environment
   \usepackage{amsmath}  % improve math presentation
   \usepackage{graphicx} % takes care of graphic including machinery
   \usepackage[margin=1in,a4paper]{geometry} % decreases margins
   \usepackage{cite} % takes care of citations
   \usepackage[final]{hyperref} % adds hyper links inside the generated pdf file
   \hypersetup{
   colorlinks=true,       % false: boxed links; true: colored links
   linkcolor=blue,        % color of internal links
   citecolor=blue,        % color of links to bibliography
   filecolor=magenta,     % color of file links
   urlcolor=blue         
   }
   %++++++++++++++++++++++++++++++++++++++++

