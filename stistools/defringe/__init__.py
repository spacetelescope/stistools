from .prepspec import prepspec
from .normspflat import normspflat
from .mkfringeflat import mkfringeflat
from .defringe import defringe

__doc__ = """
HST/STIS CCD Defringing Tools

Routines
--------

- `prepspec`     -- Calibrate STIS CCD G750L or G750M spectrum before defringing  
- `normspflat`   -- Normalize STIS CCD fringe flat  
- `mkfringeflat` -- Match fringes in STIS fringe flat to those in science data  
- `defringe`     -- Defringe by dividing the science spectrum by the fringe flat  

.. note::
    See the stistools documentation on readthedocs for the latest information and user guides:
    <https://stistools.readthedocs.io/en/latest/>
    See section 3.5.5 of the STIS Data Handbook (DHB) for more details on the defringing process:
    <http://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/stis/documentation/_documents/stis_dhb.pdf>

.. warning::
    These routines are based on PyRAF stsdas.hst_calib.stis defringing tasks, though 
    users should expect numerical discrepancies between these two implementations.
"""
