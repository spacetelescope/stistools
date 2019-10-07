#! /usr/bin/env python

import os
import re
import numpy as np
import argparse
import textwrap
from astropy.io import fits
import stistools

stistools.r_util.NOT_APPLICABLE = 'n/a'  # Fix a bug in expandFileName()

from stistools.r_util import expandFileName
# from ..basic2d import basic2d
# from ..x2d import x2d
# from ..ocrreject import ocrreject


def prepspec(inspec, outroot, darkfile=None, pixelflat=None, initguess='min'):
    """Correct STIS CCD G750L or G750M spectrum for fringing.

    Based on PyRAF `stsdas.hst_calib.stis.prepspec task 
    <https://github.com/spacetelescope/stsdas/blob/master/stsdas/pkg/hst_calib/stis/prepspec.cl>`_.

    Parameters
    ----------
    inspec: str
        Name of input 'raw' science spectrum
    outroot: str
        Root for output file name
    darkfile: str or None
        Name of superdark image.  If None, use DARKFILE in main header of input spectrum.
    pixelflat: str or None
        Name of pixel-to-pixel flat.  If None, use PIXELFLAT in main header of input
        spectrum.
    initguess: str {'min', 'med'}
        Method for initial value estimate for ocrreject.  (Default='min')

    Returns
    -------
    outname: str
        Fully qualified name of prepared spectrum (CRJ or SX2 file)

    """
    # These notes are based on the old STSDAS algorithm.  It may make sense to change the 
    # order or processing when converting to Python.

    # Note that stistools.basic2d.basic2d(), stistools.ocrreject.ocrreject(), and 
    # stistools.x2d.x2d() likely ignore ext=0 calibration flags.

    # Check inputs:
    #    - inspec:  HST/STIS, G750M/G750L, correct filetype, not already corrected
    #    - outroot:  handle dir vs filename; ".fits" in name or not
    #    - don't allow inspec to be of type {crj, sx2}
    #    - print/log name of resolved reference files
    #    - read in reference file data

    # Check if we can we CR-combine the data:
    #    - Check inspec has more than 1 {SCI, ERR, DQ} ext groups (NEXTEND//3).
    #    - Handle CR-rejection header keywords when NRPTEXP != NEXTEND//3.

    # perform only {DQICORR, BLEVCORR} via stistools.basic2d.basic2d() --> Produces temporary file
    # perform CRCORR via stistools.ocrreject.ocrreject() --> Produces CRJ file
    # perform only {BIASCORR, DARKCORR, FLATCORR, PHOTCORR; STATFLAG} via stistools.basic2d.basic2d()
    # If G750L:
    #    return CRJ file
    # If G750M:
    #    perform only {WAVECORR, HELCORR, X2DCORR} via stistools.x2d.x2d()
    #    return SX2 file
    science_data = os.path.abspath(expandFileName(inspec))
    sci_root = re.split('\.fits.*', os.path.basename(science_data), flags=re.IGNORECASE)[0].rsplit('_', 1)[0]
    opt_elem = fits.getval(science_data, 'OPT_ELEM').strip().upper()
    if (fits.getval(science_data, ext=0, keyword='DETECTOR').strip().upper() != 'CCD') or \
            (fits.getval(science_data, ext=0, keyword='INSTRUME').strip() != 'STIS'):
        raise Exception('prepspec:  Intended for use on STIS/CCD data!')


    # Make sure the necessary header keywords are set to PERFORM
    fits.setval(science_data, 'BLEVCORR', value='PERFORM')  # overscan correction
    fits.setval(science_data, 'CRCORR', value='PERFORM')  # cosmic ray rejection
    fits.setval(science_data, 'BIASCORR', value='PERFORM')  # bias subtraction
    fits.setval(science_data, 'DARKCORR', value='PERFORM')  # dark subtraction
    fits.setval(science_data, 'FLATCORR', value='PERFORM')  # I believe this is the correct flat field correction keyword

    # A few extra calibration steps if G750M data:
    if opt_elem == 'G750M':
        fits.setval(science_data, 'WAVECORR', value='PERFORM')  # wavelength calibration
        fits.setcal(science_data, 'GEOCORR', value='PERFORM')  # geometric correction
        fits.setval(science_data, 'X2DCORR', value='PERFORM')  # 2D rectification

    if darkfile is not None:
        # Use the input-specified dark file:
        fits.setval(science_data, ext=0, keyword='DARKFILE', value=dark)
    else:
        darkfile = fits.getval(science_data, ext=0, keyword='DARKFILE')
        print('Using DARKFILE from the header:\n    {}'.format(darkfile))
    # dark = os.path.abspath(expandFileName(dark))  # Expand IRAF and UNIX $VARS

    if pixelflat is not None:
        # Use the input-specified pixel flat:
        fits.setval(science_data, ext=0, keyword='PFLTFILE', value=pixelflat)
    else:
        pixelflat = fits.getval(science_data, ext=0, keyword='PFLTFILE')
        print('Using PFLTFILE from the header:\n    {}'.format(pixelflat))
    # pixflat = os.path.abspath(expandFileName(pixflat))  # Expand IRAF and UNIX $VARS

    # Calibrate!
    trl_file = '{}_trl.txt'.format(sci_root)  # Output log goes here
    cwd = os.getcwd()
    try:
        # Run calstis from within the directory with the data to find EPC files properly:
        os.chdir(os.path.dirname(science_data))
        res = stistools.calstis.calstis(os.path.basename(science_data), trailer=trl_file)
    finally:
        os.chdir(cwd)

    # Print out calstis log:
    trl_file = os.path.join(os.path.dirname(science_data), trl_file)
    with open(trl_file) as trl:
        for line in trl:
            print ('    ' + line.rstrip())

    # Raise exception on bad calstis exit code:
    if res != 0:
        raise Exception('CalSTIS exited with code {}'.format(res))

    # raise NotImplementedError()


def call_prepspec():
    """Command line entry point for prepspec().
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Correct STIS CCD G750L or G750M spectrum for fringing')
    parser.add_argument('inspec', type=str, help='Name of input "raw" science spectrum')
    parser.add_argument('outroot', type=str, help='Root for output file name')
    parser.add_argument('--darkfile', '-d', type=str, default=None,
        help='Name of superdark image. If omitted, use DARKFILE in main header of input '
             'spectrum.')
    parser.add_argument('--pixelflat', '-f', type=str, default=None,
        help='Name of pixel-to-pixel flat.  If omitted, use PIXELFLAT in main header of '
             'input spectrum.')
    parser.add_argument('--initguess', '-i', type=str, default='min', choices=['min', 'med'],
        help='Method for initial value estimate for ocrreject (Default="min")')

    args = vars(parser.parse_args())
    prepspec(**args)


if __name__ == '__main__':
    call_prepspec()
