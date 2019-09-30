#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits

from ..r_util import expandFileName


def normspflat(inflat, outflat, do_cal=True, biasfile=None, darkfile=None, pixelflat=None, 
               wavecal=None):
    """Normalize STIS CCD fringe flat.

    Based on PyRAF `stsdas.hst_calib.stis.normspflat task 
    <https://github.com/spacetelescope/stsdas/blob/master/stsdas/pkg/hst_calib/stis/normspflat.cl>`_.

    Parameters
    ----------
    inflat: str
        Name of input fringe flat
    outflat: str
        Name of normalized fringe flat output
    do_cal: bool
        Perform bias and dark subtraction and CR rejection?  (Default=True)
    biasfile: str or None
        Name of superbias image.  If None, use BIASFILE in main header of the inflat.
    darkfile: str or None
        Name of superdark image.  If None, use DARKFILE in main header of the inflat.
    pixelflat: str or None
        Name of pixel-to-pixel flat.  If None, use PFLTFILE in main header of the inflat.
    wavecal: str or None
        Name of wavecal file [ONLY FOR G750M SPECTRA].  If None, use WAVECAL in main 
        header of the inflat.

    Returns
    -------
    outname: str
        Fully qualified name of the outflat

    """
    # These notes are based on the old STSDAS algorithm.

    # strip extension from rootname of inflat (assumed to be {RAW, CRJ, SX2, X2D, clff})

    # identify wavecal

    # handle outflat filename formatting (e.g. a directory is given)

    # if do_cal:
    #    identify biasfile, darkfile, pixelflat
    #    confirm inflat has >= 2 SCI exts to allow CR-rejection
    #    handle NRPTEXP keyword
    #    PERFORM:  DQICORR, BLEVCORR, BIASCORR, DARKCORR, FLATCORR, CRCORR
    #    Update BIASFILE, DARKFILE, PFLTFILE
    #    if G750M:
    #        PERFORM:  WAVECORR, HELCORR, X2DCORR
    #    elif G750L:
    #        OMIT:  WAVECORR, HELCORR, X2DCORR
    #    call stistools.calstis.calstis(inflat, wavecal, ...)
    # 
    #    perform pixel-to-pixel flat-fielding:
    #        if G750M:
    #            make unity file from inflat's SX2 (trivial flat)
    #        elif G750L:
    #            make unity file from inflat's CRJ --> tmp
    #            copy CRJ to clff

    # elif not do_cal:
    #    if G750M:
    #        make unity file from inflat's RAW --> tmp
    #    elif G750L:
    #        make unity file from inflat's RAW --> tmp
    #        copy inflat's RAW --> clff

    # Check that G750M inflats are either sx2 or x2d.

    # Do a line-by-line cubic spline fit to the fringe flat to remove the lamp function...

    # If short-slit aperture (does not start with "52X"):
    #    Find the row with max counts in the short-slit fringe flat
    #    Uses central 60% of columns
    #    Does some flux-filtering

    # Set rows (startrow, lastrow) to be fit according to aperture name (and possibly OPT_ELEM):
    # '0.3X0.09', '0.2X0.06', '52X...'

    # Details of spline fit determined according to OPT_ELEM + CENWAVE.
    # G750M (various CENWAVEs), G750L (i.e. CENWAVE == 7751; various binning), other (never used?)
    # If G750L:
    #    Iterate spline fit

    # Put spline results into tmp file and rename to output file (passes along proper header contents)

    raise NotImplementedError()


def call_normspflat():
    """Command line entry point for normspflat().
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Normalize STIS CCD fringe flat')
    parser.add_argument('inflat', type=str, help='Name of input fringe flat')
    parser.add_argument('outflat', type=str, help='Name of normalized fringe flat output')
    parser.add_argument('--skip_cal', '-s', dest='do_cal', action='store_false', 
        help='Skip bias and dark subtraction and CR rejection?')
    parser.add_argument('--biasfile', '-b', type=str,
        help='Name of superbias image. If omitted, use BIASFILE in main header of the '
             'inflat.')
    parser.add_argument('--darkfile', '-d', type=str, 
        help='Name of superdark image. If omitted, use DARKFILE in main header of the '
             'inflat.')
    parser.add_argument('--pixelflat', '-p', type=str, 
        help='Name of pixel-to-pixel flat.  If omitted, use PFLTFILE in main header of '
             'the inflat.')
    parser.add_argument('--wavecal', '-w', type=str, 
        help='Name of wavecal file [ONLY FOR G750M SPECTRA]. If omitted, use WAVECAL in '
             'main header of the inflat.')

    args = vars(parser.parse_args())
    normspflat(**args)


if __name__ == '__main__':
    call_normspflat()
