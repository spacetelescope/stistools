#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits

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

    

    raise NotImplementedError()


def call_normspflat():
    ''' Command line entry point for normspflat().
    '''
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
