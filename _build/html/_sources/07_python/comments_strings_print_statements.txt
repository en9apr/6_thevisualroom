========================================
 Comments, Strings and Print Statements
========================================

Comments
--------

-  Non-computational parts of the program that textually describe the
   behaviour of the program.
-  Comments begin with ``#`` everything to right of the hash is ignored
   by Python.
-  Comments should be frequent so you and others can understand the
   code.

.. testcode::

    # Did you notice that these words are green? Any words after
    # a Hash Tag (#) on the same line are called comments.
    
    # Comments are lines of code that the computer does not
    # look at. These can be used to clarify the purpose of your
    # code and make it easier to read.

Strings
-------

-  Sequence of characters enclosed by a pair of single or double quotes
-  Examples are "cats hate dogs" and 'Strings are fun!'.
-  Strings are one kind of data in Python. Their data type is denoted
   ``str``.

.. testcode::

   print "This is a string"

.. testoutput::

    This is a string

.. testcode::

   print 'This is also a string'

.. testoutput::

    This is also a string

Print Statements
----------------

To print these strings and other things to the box on the right, the
word print is used.

.. testcode::

    print "Output Number One"
    print 5.069393

.. testoutput::

    Output Number One
    5.069393


Blank Lines
~~~~~~~~~~~

Print statements on their own can be used to print blank lines to the
screen, which separates output and makes it easier to read.

.. testcode::

   print 
   print "Hello" 
   print

.. testcode::
  
   <BLANKLINE>
   Hello
   <BLANKLINE>

Strings on the Same Line
~~~~~~~~~~~~~~~~~~~~~~~~

Multiple strings can be printed on the same line by using commas to
separate them. This automatically inserts a space in between the two
strings.

.. testcode::

    print "One", "Two"
    print "One","Two","Three"

.. testoutput::

    One Two
    One Two Three
    
Strings within Strings
~~~~~~~~~~~~~~~~~~~~~~

If you want to include a quotation mark or apostrophe in your string,
you need to make the symbols around the string be the opposite type.

.. testcode::

    print "'You're awesome!' he said"
    print '"Thank you!" I replied.'

.. testoutput::

    'You're awesome!' he said
    "Thank you!" I replied.


Errors with String Value
~~~~~~~~~~~~~~~~~~~~~~~~

-  If you accidentally use the same ones, errors can occur.
-  Notice that the words are not coloured in the following examples;
   that is a big clue that something is wrong.

.. testcode::

    print 'It's mine'
    print "I said "hi" to him"

.. testoutput::

      File "<ipython-input-15-a5c49eea454d>", line 1
        print 'It's mine'
                  ^
    SyntaxError: invalid syntax

Errors with Print Command
~~~~~~~~~~~~~~~~~~~~~~~~~

You can also get an error by misspelling 'print.' Notice that again,
print shows up as a different colour than normal.

.. testcode::

    prit "Error"
    printer "Error"

.. testoutput::

      File "<ipython-input-16-abc799b605cb>", line 1
        prit "Error"
                   ^
    SyntaxError: invalid syntax

Concatenate Multiple Lines on Print Command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can concatenate the output of a print command over multiple lines of
code using the comma `,` and on the same line using the plus sign `+`

.. testcode::

    print "Something is " + str(42) + " and",    # Notice that the comma is used
    print "Something else is " + str(60)

.. testoutput::

    Something is 42 and Something else is 60

-  Note: ``str`` does not add a space before and after the string
-  Commas within and between print statements do add a space

.. testcode::

    print "something",
    print "something else"
    print "--"
    print "this","and","that"

.. testoutput::

    something something else
    --
    this and that


