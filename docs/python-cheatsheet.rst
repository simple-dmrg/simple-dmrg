=================
Python cheatsheet
=================

[designed specifically for understanding and modifying simple-dmrg]

For a programmer, the standard, online `Python tutorial
<http://docs.python.org/3/tutorial/>`_ is quite nice.  Below, we try
to mention a few things so that you can get acquainted with the
``simple-dmrg`` code as quickly as possible.

Python includes a few powerful internal data structures (lists,
tuples, and dictionaries), and we use ``numpy`` (numeric python) and
``scipy`` (additional "scientific" python routines) for linear
algebra.

Basics
------

Unlike many languages where blocks are denoted by braces or special
``end`` statements, blocks in python are denoted by indentation level.
Thus indentation and whitespace are significant in a python program.

It is possible to execute python directly from the commandline::

    $ python

This will bring you into python's real-eval-print loop (REPL).  From
here, you can experiment with various commands and expressions.  The
examples below are taken from the REPL, and include the prompts
("``>>>``" and "``...``") one would see there.

Lists, tuples, and loops
------------------------

The basic sequence data types in python are lists and tuples.

A ``list`` can be constructed literally::

    >>> x_list = [2, 3, 5, 7]

and a number of operations can be performed on it::

    >>> len(x_list)
    4

    >>> x_list.append(11)
    >>> x_list
    [2, 3, 5, 7, 11]

    >>> x_list[0]
    2

    >>> x_list[0] = 0
    >>> x_list
    [0, 3, 5, 7, 11]

Note, in particular, that python uses indices counting from zero, like C (but unlike Fortran and Matlab).

A ``tuple`` in python acts very similarly to a list, but once it is constructed it cannot be modified.  It is constructed using parentheses instead of brackets::

    >>> x_tuple = (2, 3, 5, 7)

Lists and tuples can contain any data type, and the data type of the elements need not be consistent::

    >>> x = ["hello", 4, 8, (23, 12)]

It is also possible to get a subset of a list (e.g. the first three
elements) by using Python's `slice notation
<http://stackoverflow.com/questions/509211/pythons-slice-notation>`_::

    >>> x = [2, 3, 5, 7, 11]
    >>> x[:3]
    [2, 3, 5]

Looping over lists and tuples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Looping over a ``list`` or ``tuple`` is quite straightforward::

    >>> x_list = [5, 7, 9, 11]
    >>> for x in x_list:
    ...     print(x)
    ... 
    5
    7
    9
    11

If you wish to have the corresponding indices for each element of the
list, the ``enumerate()`` function will provide this::

    >>> x_list = [5, 7, 9, 11]
    >>> for i, x in enumerate(x_list):
    ...     print(i, x)
    ... 
    0 5
    1 7
    2 9
    3 11

If you have two (or more) parallel arrays with the same number of
elements and you want to loop over each of them at once, use the
``zip()`` function::

    >>> x_list = [2, 3, 5, 7]
    >>> y_list = [12, 13, 14, 15]
    >>> for x, y in zip(x_list, y_list):
    ...     print(x, y)
    ... 
    2 12
    3 13
    5 14
    7 15

There is a syntactic shortcut for transforming a list into a new one,
known as a `list comprehension <http://docs.python.org/3/tutorial/datastructures.html#list-comprehensions>`_::

    >>> primes = [2, 3, 5, 7]
    >>> doubled_primes = [2 * x for x in primes]
    >>> doubled_primes
    [4, 6, 10, 14]

Dictionaries
------------

Dictionaries are python's powerful mapping data type.  A number,
string, or even a tuple can be a key, and any data type can be the
corresponding value.

Literal construction syntax::

    >>> d = {2: "two", 3: "three"}

Lookup syntax::

    >>> d[2]
    'two'
    >>> d[3]
    'three'

Modifying (or creating) elements::

    >>> d[4] = "four"
    >>> d
    {2: 'two', 3: 'three', 4: 'four'}

The method ``get()`` is another way to lookup an element, but returns
the special value ``None`` if the key does not exist (instead of
raising an error)::

    >>> d.get(2)
    'two'
    >>> d.get(4)

Looping over dictionaries
~~~~~~~~~~~~~~~~~~~~~~~~~

Looping over the keys of a dictionary::

    >>> d = {2: "two", 3: "three"}
    >>> for key in d:
    ...     print(key)
    ... 
    2
    3

Looping over the values of a dictionary::

    >>> d = {2: "two", 3: "three"}
    >>> for value in d.values():
    ...     print(value)
    ... 
    two
    three

Looping over the keys and values, together::

    >>> d = {2: "two", 3: "three"}
    >>> for key, value in d.items():
    ...     print(key, value)
    ... 
    2 two
    3 three

Functions
---------

Function definition in python uses the ``def`` keyword::

    >>> def f(x):
    ...     y = x + 2
    ...     return 2 * y + x
    ... 

Function calling uses parentheses, along with any arguments to be passed::

    >>> f(2)
    10
    >>> f(3)
    13

When calling a function, it is also possibly to specify the arguments by name (e.g. ``x=4``)::

    >>> f(x=4)
    16

An alternative syntax for writing a one-line function is to use python's ``lambda`` keyword::

    >>> g = lambda x: 3 * x
    >>> g(5)
    15

numpy arrays
------------

``numpy`` provides a multi-dimensional array type.  Unlike lists and
tuples, ``numpy`` arrays have fixed size and hold values of a single
data type.  This allows the program to perform operations on large
arrays very quickly.

Literal construction of a 2x2 matrix::

    >>> np.array([[1, 2], [3, 4]], dtype='d')
    array([[ 1.,  2.],
	   [ 3.,  4.]])

Note that ``dtype='d'`` specifies that the type of the array should
be double-precision (real) floating point.

It is also possibly to construct an array of all zeros::

    >>> np.zeros([3, 4], dtype='d')
    array([[ 0.,  0.,  0.,  0.],
	   [ 0.,  0.,  0.,  0.],
	   [ 0.,  0.,  0.,  0.]])

And then elements can be added one-by-one::

    >>> x = np.zeros([3, 4], dtype='d')
    >>> x[1, 2] = 12
    >>> x[1, 3] = 18
    >>> x
    array([[  0.,   0.,   0.,   0.],
	   [  0.,   0.,  12.,  18.],
	   [  0.,   0.,   0.,   0.]])

It is possible to access a given row or column by index::

    >>> x[1, :]
    array([  0.,   0.,  12.,  18.])
    >>> x[:, 2]
    array([  0.,  12.,   0.])

or to access multiple columns (or rows) at once::

    >>> col_indices = [2, 1, 3]
    >>> x[:, col_indices]
    array([[  0.,   0.,   0.],
	   [ 12.,   0.,  18.],
	   [  0.,   0.,   0.]])

For matrix-vector (or matrix-matrix) multiplication use the
``np.dot()`` function::

    >>> np.dot(m, v)

.. warning::

    One tricky thing about ``numpy`` arrays is that they do not act as
    matrices by default.  In fact, if you multiply two ``numpy``
    arrays, python will attempt to multiply them element-wise!

To take an inner product, you will need to take the
transpose-conjugate of the left vector yourself::

    >>> np.dot(v1.conjugate().transpose(), v2)

Array storage order
~~~~~~~~~~~~~~~~~~~

Although a ``numpy`` array acts as a multi-dimensional object, it is
actually stored in memory as a one-dimensional contiguous array.
Roughly speaking, the elements can either be stored column-by-column
("column major", or "Fortran-style") or row-by-row ("row major", or
"C-style").  As long as we understand the underlying storage order of
an array, we can reshape it to have different dimensions.  In
particular, the logic for taking a partial trace in ``simple-dmrg``
uses this reshaping to make the system and environment basis elements
correspond to the rows and columns of the matrix, respectively.  Then,
only a simple matrix multiplication is required to find the reduced
density matrix.

Mathematical constants
----------------------

``numpy`` also provides a variety of mathematical constants::

    >>> np.pi
    3.141592653589793
    >>> np.e
    2.718281828459045

Experimentation and getting help
--------------------------------

As mentioned above, python's REPL can be quite useful for
experimentation and getting familiar with the language.  Another thing
we can do is to import the ``simple-dmrg`` code directly into the REPL
so that we can experiment with it directly.  The line::

    >>> from simple_dmrg_01_infinite_system import *

will execute all lines *except* the ones within the block that says::

    if __name__ == "__main__":

So if we want to use the finite system algorithm, we can (assuming our
source tree is in the ``PYTHONPATH``, which should typically include
the current directory)::

    $ python
    >>> from simple_dmrg_04_eigenstate_prediction import *
    >>> finite_system_algorithm(L=10, m_warmup=8, m_sweep_list=[8, 8, 8])

It is also possible to get help in the REPL by using python's built-in
``help()`` function on various objects, functions, and types::

    >>> help(sum)   # help on python's sum function

    >>> help([])    # python list methods
    >>> help({})    # python dict methods

    >>> help({}.setdefault)   # help on a specific dict method

    >>> import numpy as np
    >>> help(np.log)          # natural logarithm
    >>> help(np.linalg.eigh)  # eigensolver for hermitian matrices
