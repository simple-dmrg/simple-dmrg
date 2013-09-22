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

Once the relevant software is installed, each program is contained
entirely in a single file.  The first program, for instance, can be
run by issuing::

    $ python simple_dmrg_01_infinite_system.py

.. note::

    If you see an error that looks like this::

        SyntaxError: future feature print_function is not defined

    then you are using a version of Python below 2.6.  Although it
    would be best to upgrade, it may be possible to make the code work
    on Python versions below 2.6 without much trouble.
