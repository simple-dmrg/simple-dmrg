===========
SIMPLE DMRG
===========

Source code: https://github.com/simple-dmrg/simple-dmrg/

Documentation: http://simple-dmrg.readthedocs.org/

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
the included :doc:`Python cheatsheet <python-cheatsheet>`, which lists
many of the Python features used by ``simple-dmrg``, and which should
be helpful when trying the included :doc:`exercises <exercises>`.)

The four modules build up DMRG from its simplest implementation to
more complex implementations and optimizations.  Each file adds lines
of code and complexity compared with the previous version.

1. :doc:`Infinite system algorithm <01_infinite_system>`
   (~180 lines, including comments)
2. :doc:`Finite system algorithm <02_finite_system>`
   (~240 lines)
3. :doc:`Conserved quantum numbers <03_conserved_quantum_numbers>`
   (~310 lines)
4. :doc:`Eigenstate prediction <04_eigenstate_prediction>`
   (~370 lines)

Throughout the tutorial, we focus on the spin-1/2 `Heisenberg XXZ
model <http://en.wikipedia.org/wiki/Heisenberg_model_(quantum)>`_, but
the code could easily be modified (or expanded) to work with other
models.

Authors
=======

- James R. Garrison (UCSB)
- Ryan V. Mishmash (UCSB)

Licensed under the MIT license.  If you plan to publish work based on
this code, please contact us to find out how to cite us.

Contents
========

.. toctree::
   :maxdepth: 2

   using
   exercises
   python-cheatsheet
   references
   source-code
