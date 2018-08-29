.. stistools documentation master file, created by
   Warren Hack on Mon Oct 1 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========================================
Welcome to stistool's documentation!
========================================
Stistools is a package that provides Python-based data processing tools for working with Space Telescope Imaging Spectrograph
(STIS) data. It contains the full STIS calibration pipeline as well as its individual components should the user wish to do their
calibrations manually. Additionally, stistools features a selection of analysis tools independent from the pipeline. These analysis
tools are in active development, so there is more to come.

Relationship to the IRAF/PyRAF stistools Package
*******************************************

The intent of this package is to provide the user community with many of the useful tools contained within the original IRAF/PyRAF
stistools package in Python. This is consistent with a larger effort to transition the user community from IRAF/PyRAF to Python,
which can be read about `here <https://stak-notebooks.readthedocs.io/en/latest/index.html>`_. At this time, there remains several tools
from the original IRAF/PyRAF version of the package that are still in development in the Python package. In this transitional period,
it’s recommended that STIS users familiarize themselves with this package. Understanding that many of the users of this package
will be those who may be newer to Python as a result of this transition, this documentation seeks to provide a step-by-step guide
for getting started.

========================
Getting Started
========================
.. toctree::
   :maxdepth: 2

   Installation and Setup <gettingstarted.rst>
=======================
Pipeline Software
=======================
.. toctree::
   :maxdepth: 2

   basic2d
   x1d
   calstis
   wavelen
   ocrreject

=======================
Analysis Tools
=======================
.. toctree::
   :maxdepth: 2

   sshift <sshift.rst>
   stisnoise <stisnoise.rst>
   mktrace <mktrace.rst>
   evaldisp <evaldisp.rst>
   wx2d
   radialvel
   r_util
   gettable
   inttag

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
