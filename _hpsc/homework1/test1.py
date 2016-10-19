"""
test1.py file

Note: this is a Python script.  What you are reading now in the triple
quotes is a comment (a "docstring" since it's the first thing in this file).

This file is to test whether you have environment variables set properly and
are able to import numpy (for numerical computation), matplotlib, 
and pylab (for plotting).

Instructions:

Make sure your environment variables $UWHPSC and $MYHPSC are set properly.

Insert your name where instructed below.

Type
    $ python test1.py
at the Unix prompt ($ represents the prompt) to run the code.

If the output looks ok, change the file to set debug_mode = False
below and run it again.  This time the results will go into a file
test1output.txt.   Check that these look ok.

Commit the final version of this file along with test1output.txt to your git
repository.

"""

import os,sys

debug_mode = False

if debug_mode:
    output_file = sys.stdout  # send output to screen
else:
    # Once it's working, send output to a file rather than screen:
    output_file = open("test1output.txt","w")
    sys.stdout = output_file

print "Code run by **Andrew Roberts**"

UWHPSC = os.environ.get('UWHPSC','** not set **')
print "Environment variable UWHPSC is %s"  % UWHPSC

MYHPSC = os.environ.get('MYHPSC','** not set **')
print "Environment variable MYHPSC is %s"  % MYHPSC


# Check for various Python modules and print an error message if
# an exception occurs:

try:
    import numpy
    print "Imported numpy ok"
except:
    print "** Error importing numpy"

try:
    import matplotlib
    print "Imported matplotlib ok"
except:
    print "** Error importing matplotlib"

try:
    import pylab
    print "Imported pylab ok"
except:
    print "** Error importing pylab"

if not debug_mode:
    output_file.close()  # close the text file


