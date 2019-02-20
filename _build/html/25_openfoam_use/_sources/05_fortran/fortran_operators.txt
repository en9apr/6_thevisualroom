=========
Operators
=========

.. contents::
   :local:

.. highlight:: fortran

Associativity
=============

* Right-to-left means the value on the right is evaluated first
* Left-to-right means the value on the left is evaluated first

Order of Precedence
===================

The order of precedence for operators:

* Bracket
* Arithmetic
* Relational
* Logical
* Assignment

Bracket Operators
-----------------

.. list-table::
   :header-rows: 1
   :widths: 10 10

   * - Operator
     - Meaning
   * - ::

           (

     - Left bracket
   * - ::

           )

     - Right bracket

Arithmetic Operators
--------------------

All operators here are binary (require two numbers):

Only ``+`` and ``-`` are both binary and unary operators.

.. list-table::
   :header-rows: 1
   :widths: 10 10 15 15

   * - Precedence
     - Operator
     - Meaning
     - Associativity
   * - 1
     - ::

           **

     - Exponentiation
     - Right-to-left
   * - =2
     - ::

           *
     - Muliplication
     - Left-to-right
   * - =2
     - ::

           /
     - Division
     - Left-to-right
   * - =3
     - ::

           +
     - Addition
     - Left-to-right
   * - =3
     - ::
 
           -
     - Subtraction
     - Left-to-right


Relational Operators
--------------------

All are binary. As well as numbers, you can also compare strings, e.g. ``"string_1" .ne. "string_2"`` or ``"string_1" .eq. "string_2"``

.. list-table::
   :header-rows: 1
   :widths: 10 10 10 15 15

   * - Precedence
     - Operator 1
     - Operator 2
     - Meaning
     - Associativity
   * - =1
     - ::

           .lt.
     - ::

           <
     - Less Than
     - None
   * - =1
     - ::

           .le.
     - ::

           <=
     - Less Than or Equal To
     - None
   * - =1
     - ::

           .eq.
     - ::

           ==
     - Equal To
     - None
   * - =1
     - ::

           .ne.
     - ::

           /=
     - Not Equal To
     - None
   * - =1
     - ::
 
           .gt.
     - ::

           >
     - Greater Than
     - None
   * - =1
     - ::

           .ge.
     - ::

           >=
     - Greater Than or Equal To
     - None

Logical Operators
-----------------

Logical Not is unary, the rest are binary:

.. list-table::
   :header-rows: 1
   :widths: 10 10 15 15

   * - Precedence
     - Operator
     - Meaning
     - Associativity
   * - 1
     - ::

           .not.

     - Unary Logical Not
     - Right-to-left
   * - 2
     - ::

           .and.
     - Binary Logical And
     - Left-to-right
   * - 3
     - ::

           .or.
     - Binary Logical Or
     - Left-to-right
   * - =4
     - ::

           .eqv.
     - Binary Logical Equivalence
     - Left-to-right
   * - =4
     - ::
 
           .neqv.
     - Binary Logical Not Equivalence
     - Left-to-right

Assignment Operators
--------------------

Assignment operator is binary:

.. list-table::
   :header-rows: 1
   :widths: 10 10 15

   * - Operator
     - Meaning
     - Associativity
   * - ::

           =

     - Assignment
     - Right-to-left

Concatenation, Continuation and Comments
========================================

.. list-table::
   :header-rows: 1
   :widths: 10 15

   * - Operation
     - Meaning
   * - ::

           "text_1" // "text_2"

     - String concatenation
   * - ::

           A = 175.5 * Year &
               + Count / 100
     - Ampersand (&) is a continuation line
   * - ::

           ! Text
     - Exclamation Mark (!) is a comment
