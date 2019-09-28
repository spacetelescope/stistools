from .prepspec import prepspec
from .normspflat import normspflat
from .mkfringeflat import mkfringeflat
from .defringe import defringe

__doc__ = """
HST/STIS CCD Defringing Tools

Routines
--------

- `prepspec`     -- Correct STIS CCD G750L or G750M spectrum for fringing  
- `normspflat`   -- Normalize STIS CCD fringe flat  
- `mkfringeflat` -- Match fringes in STIS fringe flat to those in science data  
- `defringe`     -- Correct STIS CCD G750L or G750M spectrum for fringing  

.. note::
    See section 3.5.5 of the `STIS Data Handbook (DHB) <http://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/stis/documentation/_documents/stis_dhb.pdf>`_
    for more details on the defringing process.

.. warning::
    These routines are based on PyRAF stsdas.hst_calib.stis defringing tasks, though 
    users should expect numerical discrepancies between these two implementations.

Example
-------

.. code-block:: python

    from stistools import defringe
    from astropy.io import fits
    
    # Apply preliminary calibrations to the raw science data (if necessary):
    # (The output CRJ may already exist in pipeline products.)
    defringe.prepspec('o8u201060_raw.fits', 'o8u201060_crj.fits')  # G750L
    
    # Normalize the contemporaneous fringe flat:
    defringe.normspflat('o8u201070_raw.fits', 'o8u201070_nsp.fits, do_cal=True)
    
    # Remove fringes due to order-sorter filter, as these are already included in the
    # sensitivity function:
    with fits.open('o8u201070_nsp.fits', 'update') as f:
        f[1].data[:,:250] = 1.
    
    # Fit the normalized fringe flat to the science data in shift and amplitude:
    defringe.mkfringeflat('o8u201060_crj.fits', 'o8u201070_nsp.fits', 'o8u201070_frr.fits')
    
    # Apply the fitted fringe flat to the 
    defringe.defringe('o8u201060_crj.fits', 'o8u201070_frr.fits', 'o8u201060_drj.fits')
"""
