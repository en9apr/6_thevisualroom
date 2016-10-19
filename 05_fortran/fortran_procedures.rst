====================
Intrinsic Procedures
====================

.. contents::
   :local:

.. highlight:: fortran

Transfer and Conversion Functions
=================================

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Code
     - Meaning
   * - ::

             abs(a) 
     - Absolute value
   * - ::

             aimag(z)
     - Imaginary part of complex z
   * - ::

             aint(x, kind), anint(x, kind)
     - To whole number real
   * - ::

             dble(a) 
     - To double precision
   * - ::

             cmplx(x,y, kind)
     - Create x + iy (y optional)
   * - ::
          
             int(a, kind), nint(a, kind)
     - To int (truncated/rounded)
   * - ::
 
             real(x, kind)
     - To real
   * - ::

             conj(z)
     - Complex conjugate
   * - ::

             char(i, kind), achar(i)
     - Char of ASCII code (pure 7bit)
   * - ::

             ichar(c), iachar(c)
     - ASCII code of character
   * - ::
 
             logical(l, kind)
     - Change kind of logical l 
   * - ::

             ibits(i, pos, len)
     - Extract sequence of bits
   * - ::

             transfer(source, mold, size)
     - Reinterpret data


Arrays and Matrices
===================

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Code
     - Meaning
   * - ::

             allocated(a)
     - Check if array is allocated
   * - ::
   
             lbound(a, dim), ubound(a,dim)
     - Lowest/highest index in array
   * - ::

             shape(a)
     - Shape (dimensions) of array
   * - ::

             size(array, dim)
     - Extent of array along dim
   * - ::

             all(mask, dim), any(mask, dim)
     - Check boolean array
   * - ::

             count(mask, dim)
     - Number of true elements
   * - ::

             maxval(a,d,m), minval(a,d,m)
     - Find max/min in masked array
   * - ::

             product(a, dim, mask)
     - Product along masked dimen.
   * - ::

             dot_product(a, b)
     - Dot product
   * - ::

             sum(a, dim, mask)
     - Sum along masked dimension
   * - ::

             merge(tsource, fsource, mask)
     - Combine arrays as mask says   
   * - ::

             pack(array, mask, vector)
     - Packs masked array into vect.
   * - ::

             unpack(vector, mask , field)
     - Unpack vect. into masked field
   * - ::

             spread(source, dim, n)
     - Extend source array into dim.
   * - ::

             reshape(src,shape,pad,order)
     - Make array of shape from src
   * - ::

             cshift(a,s,d),eoshift(a,s,b,d)
     - (Circular) shift   
   * - ::

             transpose(matrix)
     - Transpose a matrix
   * - ::

             maxloc(a,mask), minloc(a,mask)
     - Find pos. of max/min in array

Computation Functions
=====================

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Code
     - Meaning
   * - ::

             ceiling(a), floor(a)
     - To next higher/lower int
   * - ::

             conj(z)
     - complex conjugate
   * - ::

             dim(x,y)
     - max(x-y, 0)
   * - ::

             max(a1, a2, a3), 
             min(a1, a2, a3)
     - maximum/minimum
   * - ::

             dprod(a,b)
     - dp product of  sp a, b
   * - ::

             mod(a,p)
     - a mod p
   * - ::

             modulo(a,p)
     - modulo with sign of a/p
   * - ::

             sign(a,b)
     - make sign of a = sign of b
   * - ::

             matmul(m1, m2)
     - matrix multiplication
   * - ::

             dot_product(a,b)
     - dot product of vectors
   * - ::

             sin(a), cos(a), tan(a),
             sinh(a), cosh(a), tanh(a)
     - Trigonometric functions where :math:`a` is in radians
   * - ::

             dsin(a), dcos(a), dtan(a), 
             dsinh(a), dcosh(a), dtanh(a)
     - Trigonometric functions where :math:`a` is in degrees
   * - ::

             acos(b), asin(b), atan(b)
     - Trigonometric functions where the result is in radians
   * - ::

             dacos(b), dasin(b), datan(b)
     - Trigonometric functions where the result is in degrees
   * - ::
             
             exp(a)
     - Exponential function
   * - ::

             log(a)
     - Natural logarithm
   * - ::

             log10(a)
     - Base 10 logarithm
   * - ::

             sqrt(a)
     - Square root


Numeric Inquiry and Manipulation Functions
==========================================

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Code
     - Meaning
   * - ::

             kind(x)
     - Kind-parameter of variable x
   * - ::

             digits(x)
     - Significant digits in model
   * - ::

             bit_size(i)
     - Number of bits for int in model
   * - ::

             epsilon(x)
     - Small pos. number in model
   * - ::

             huge(x)
     - Largest number in model
   * - ::

             minexponent(x)
     - Smallest exponent in model
   * - ::

             maxexponent(x)
     - Largest exponent in model
   * - ::

             precision(x)
     - Decimal precision for reals in
   * - ::

             radix(x)
     - Base of the model
   * - ::

             range(x)
     - Dec. exponent range in model
   * - ::

             tiny(x)
     - Smallest positive number
   * - ::

             exponent(x)
     - Exponent part of x  in model
   * - ::

             fraction(x)
     - Fractional part of x in model
   * - ::

             nearest(x)
     - Nearest machine number
   * - ::

             rrspacing(x)
     - Reciprocal of relative spacing
   * - ::

             scale(x,i) 
     - x b**i
   * - ::

             set_exponent(x,i)
     - x b**(i-e)
   * - ::

             spacing(x)
     - Absolute spacing of model

String Functions
================

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Code
     - Meaning
   * - ::

             lge(s1,s2), lgt, lle, llt
     - String comparison
   * - ::

             adjustl(s), adjustr(s)
     - left- or right-justify string
   * - ::

             index(s, sub, from_back)
     - find substr. in string (or 0)
   * - ::

             trim(s)
             trim(adjustl(s))
     - s without trailing blanks
       s without leading or trailing blanks
       
   * - ::

             len_trim(s)
     - length of s, w/ trailing blanks
   * - ::

             scan(s, setd, from_back)
     - search for any char in set
   * - ::

             verify(s, set, from_back)
     - check for presence of set-chars
   * - ::

             len(string)
     - length of string
   * - ::

             repeat(string, n)
     - concat n copies of string


Bit Functions (on integers)
===========================

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Code
     - Meaning
   * - ::

             btest(i,pos)
     - Test bit of integer value
   * - ::

             iand(i,j), ieor(i,j), ior(i,j)
     - And, xor, or of bit in 2 integers
   * - ::

             ibclr(i,pos), ibset(i, pos)
     - Set bit of integer to 0 / 1
   * - ::

             ishft(i, sh), ishftc(i, sh, s)
     - Shift bits in i
   * - ::

             not(i)
     - Bit-reverse integer

Miscellaneous Intrinsic Subroutines
===================================

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Code
     - Meaning
   * - ::

             date_and_time(d, t, z, v)
     - put current time in d,t,z,v
   * - ::

             mvbits(f, fpos, len, t, tpos)
     - copy bits between int vars
   * - ::

             random_number(harvest)
     - fill harvest randomly
   * - ::

             random_seed(size, put, get)
     - restart/query random generator
   * - ::

             system_clock(c, cr, cm)
     - get processor clock info

















