Installation and Setup
=====================

=====================
Installing stistools
=====================


Through Astroconda (Recommended)
--------------------------------
The simplest installation method for installing stistools is to install the "Standard Software Stack" of AstroConda.
AstroConda is an STScI-maintained software distribution channel for Conda, a package and environment management system.
The standard software stack of AstroConda contains all of STScI's publicly distributed software, as well as all of the
dependencies to run them. Effectively, it takes care of everything for you.

The first step is to download Conda. There are a few different flavors of Conda, but for most cases we'd recommend installing
the Python 3 version of miniconda. For a step-by-step guide on installing Conda, consult the
`AstroConda documentation <https://astroconda.readthedocs.io/en/latest/getting_started.html#installing-conda-the-choice-is-yours>`_.

With Conda installed, let's now create an AstroConda environment. For this step, you can continue to follow the
`AstroConda documentation <https://astroconda.readthedocs.io/en/latest/getting_started.html#installing-conda-the-choice-is-yours>`_,
which contains information on the different stack choices, or you can just read on and execute the commands we'd recommend for stistools.

First, we need to configure Conda to grant it access to the AstroConda channel, which can be done by running the following command in a
BASH shell.

.. code-block:: sh

    $ conda config --add channels http://ssb.stsci.edu/astroconda
    # Writes changes to ~/.condarc

Now we're ready to install the STScI Software Stack, we'll accomplish this by setting up a fresh Conda environment.
Run the following command, you can change "astroconda" to be whatever name for the environment you wish.

.. code-block:: sh

    $ conda create -n astroconda stsci

Once the installation is complete, you can access your new environment by activating it:

.. code-block:: sh

    $ source activate astroconda

Once activated, you now have access to all of the STScI software, including stistools! If you want to deactivate an environment,
you can do so like this:

.. code-block:: sh

    $ source deactivate astroconda

Keep in mind that whenever you open a new terminal, by default your environment will not be activated (this can be changed). So be sure to activate it before
attempting to use stistools. When in your environment, you can now interact with stistools like any other Python package.

.. code-block:: sh

    $ python

    >>> import stistools
    The following tasks in the stistools package can be run with TEAL:
       basic2d      calstis     ocrreject     wavecal        x1d          x2d



=====================
Setting up CRDS (Optional)
=====================
