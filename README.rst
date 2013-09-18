===========
SIMPLE DMRG
===========

On the web: https://github.com/simple-dmrg/simple-dmrg/

The goal of this tutorial (given at the `2013 summer school on quantum
spin liquids <http://www.democritos.it/qsl2013/>`_, in Trieste, Italy)
is to present the `density-matrix renormalization group
<http://en.wikipedia.org/wiki/Density_matrix_renormalization_group>`_
(DMRG) in its traditional formulation (i.e. without using matrix
product states).  DMRG is a numerical method that allows for the
efficient simulation of quantum model Hamiltonians.  Since it is a
low-entanglement approximation, it often works quite well for
one-dimensional systems, giving results that are nearly exact.

Typical implementations of DMRG in C++ or Fortran can be tens of
thousands of lines long.  Here, we have attempted to strike a balance
between clear, simple code, and including many features and
optimizations that would exist in a production code.  One thing that
helps with this is the use of `Python <http://www.python.org/>`_.  We
have tried to write the code in a very explicit style, hoping that it
will be (mostly) understandable to somebody new to Python. (See also
the included `Python cheatsheet <python-cheatsheet.rst>`_, which lists
many of the Python features used by ``simple-dmrg``, and which should
be helpful when trying the included `exercises <exercises_day1.rst>`_.)

The four modules build up DMRG from its simplest implementation to
more complex implementations and optimizations.  Each file adds lines
of code and complexity compared with the previous version.

1. `Infinite system algorithm <simple_dmrg_01_infinite_system.py>`_
   (~180 lines, including comments)
2. `Finite system algorithm <simple_dmrg_02_finite_system.py>`_
   (~240 lines)
3. `Conserved quantum numbers <simple_dmrg_03_conserved_quantum_numbers.py>`_
   (~310 lines)
4. `Eigenstate prediction <simple_dmrg_04_eigenstate_prediction.py>`_
   (~370 lines)

Throughout the tutorial, we focus on the spin-1/2 `Heisenberg XXZ
model <http://en.wikipedia.org/wiki/Heisenberg_model_(quantum)>`_, but
the code could easily be modified (or expanded) to work with other
models.

Using the code
==============

The requirements are:

* `Python <http://www.python.org/>`_ 2.6 or higher (Python 3 works as well)
* `numpy <http://www.numpy.org/>`_ and `scipy <http://www.scipy.org/>`_

Download the code using the `Download ZIP
<https://github.com/simple-dmrg/simple-dmrg/archive/master.zip>`_
button on github, or run the following command from a terminal::

    $ wget -O simple-dmrg-master.zip https://github.com/simple-dmrg/simple-dmrg/archive/master.zip

Within a terminal, execute the following to unpack the code::

    $ unzip simple-dmrg-master.zip
    $ cd simple-dmrg-master/

.. note::

    If you are using the computers in the lab at SISSA, you will need
    to activate three modules for the correct version of Python to
    run.  These commands are given in `activate_sissa_modules
    <activate_sissa_modules>`_, or you can run the following command
    to activate them::

        $ source activate_sissa_modules

Once the relevant software is installed, each program is contained
entirely in a single file.  The first program, for instance, can be
run by issuing::

    $ python simple_dmrg_01_infinite_system.py

.. note::

    If you see an error that looks like this::

        SyntaxError: future feature print_function is not defined

    then you are using a version of Python below 2.6.  If you are in
    the SISSA computer lab, follow the instructions above to activate
    the correct modules.  If you are using your own machine, we may be
    able to make the code work on versions of Python prior to 2.6
    without much trouble.

.. note::

    The version of Python in the SISSA computer lab is compiled with
    `many debugging checks enabled
    <http://docs.python.org/2/c-api/intro.html#debugging-builds>`_.
    For instance, after running the code, you may see something that
    looks like::

        [32067 refs]

    As a result of these checks, the code runs approximately 15 times
    slower on these machines than it does on our laptops.

Authors
=======

* James R. Garrison (UCSB)
* Ryan V. Mishmash (UCSB)

Licensed under the MIT license.  If you plan to publish work based on
this code, please contact us to find out how to cite us.
