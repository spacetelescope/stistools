Installation and Setup
======================

====================================
Installing ``stistools`` (``stenv``)
====================================

The simplest installation method for installing ``stistools`` is to install the ``stenv`` environment. ``stenv`` is an
STScI-maintained software distribution for Conda, a package and environment management system. The standard software
stack of ``stenv`` contains all of STScI's publicly distributed software, as well as all of the dependencies to run
them. Effectively, it takes care of everything for you.

The first step is to download Conda. There are a few different flavors of Conda, but for most cases we'd recommend
installing the Python 3 version of ``miniconda``. For a step-by-step guide on installing Conda, consult the
`documentation <https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links>`_.

With Conda installed, let's now create an ``stenv`` environment. For this step, you can follow the
`stenv documentation <https://github.com/spacetelescope/stenv>`_. The following instructions assume you have created a
``conda`` environment called ``stisenv``.

Once the installation is complete, you can access your new environment by activating it:

.. code-block:: shell

  conda activate stisenv

Once activated, you now have access to all of the STScI software, including ``stistools``! If you want to deactivate an
environment, you can do so like this:

.. code-block:: shell

  conda deactivate stisenv

Keep in mind that whenever you open a new terminal, by default your environment will not be activated (this can be
changed). So be sure to activate it before attempting to use ``stistools``. When in your environment, you can now
interact with ``stistools`` like any other Python package.

.. code-block:: shell

  python
  >>> import stistools
  The following tasks in the stistools package can be run with TEAL:
     basic2d      calstis     ocrreject     wavecal        x1d          x2d

===========================================
Getting the Latest Version of ``stistools``
===========================================
Sometimes, it may be the case that new additions to ``stistools`` have not yet been packaged into a proper release
through ``stenv``. In these instances, the installation of ``stistools`` through ``stenv`` will not contain the most
recent additions to the package. The following instructions outline how to grab and install the latest version of
``stistools``, if you require something that has been released very recently.

To start, we'll assume that you've gone through the process above, installing ``stistools`` through ``stenv``. Even
though ``stenv`` does not contain our most up-to-date version of ``stistools`` in this case, it does still provide us
with all of the necessary dependencies needed to run ``stistools``.

First, let's clone the github ``stistools`` repository down to our local machines. This essentially downloads the latest
stable version of the package to your computer. We can clone ``stistools`` by running the following command:

.. code-block:: shell

  git clone https://github.com/spacetelescope/stistools.git

Note that this will create a "``stistools``" folder in your local directory. Navigate into this directory once the clone
finishes executing. We want to install this on top of our ``stenv`` environment, so activate your desired environment
like so:

.. code-block:: shell

  conda activate stisenv

Because developer versions of ``stistools`` share the same version numbers as the last release, we'll need to remove the
version of ``stistools`` that came with our ``stenv`` environment, we can do this through conda:

.. code-block:: shell

  conda uninstall --force stistools

The ``--force`` flag is necessary for instructing conda not to uninstall packages that depend on ``stistools``. We can
now install the latest version of ``stistools``. In the ``stistools`` directory, run:

.. code-block:: shell

  pip install .

This builds the ``stistools`` package up based on the source code we cloned to our local machines. Note that this
overwrites the existing version of ``stistools`` that was installed through ``stenv``. With this, you should now have
the latest version of ``stistools`` installed in your ``stisenv`` environment.


=============================
Setting up CRDS (Recommended)
=============================

Some calibration tasks in ``stistools`` require additional reference files to successfully run. In the past, users were
expected to download these reference files manually by using `MAST <http://archive.stsci.edu/hst/search.php>`_. While
this approach is still valid, it can be inconvenient. The HST Calibration Reference Data System (CRDS) has a
`python package <https://hst-crds.stsci.edu/docs/cmdline_bestrefs/>`_ that can easily download and cache the relevant
reference files for your data for you. And in fact, the crds package is a part of the ``stenv`` stack and therefore is
already installed if you've installed ``stistools`` through ``stenv``. To get this setup, all we need to do is run a few
commands:

.. code-block:: shell

  export CRDS_PATH="$HOME/crds_cache"
  export CRDS_SERVER_URL="https://hst-crds.stsci.edu"
  export oref="${CRDS_PATH}/references/hst/oref/"

The above syntax define where your personal copies of CRDS reference files will be stored and the CRDS server that is
used. Then the following command may be used to assign and obtain the best references files:

.. code-block:: shell

  crds bestrefs --update-bestrefs --sync-references=1 --files *.fits

Note that in this example bestrefs will run on files currently in your working directly. You can modify where it looks
by updating the final input.
