#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits

from ..r_util import expandFileName
from ..basic2d import basic2d
from ..x2d import x2d
from ..ocrreject import ocrreject


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

    raise NotImplementedError()


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
