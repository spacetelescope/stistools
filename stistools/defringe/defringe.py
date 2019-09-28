#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits

def defringe(inspec, inflat, outspec, do_cal=False, darkfile=None, pixelflat=None,
             do_defringe=True):
    """Correct STIS CCD G750L or G750M spectrum for fringing.

    Based on PyRAF `stsdas.hst_calib.stis.defringe task
    <https://github.com/spacetelescope/stsdas/blob/master/stsdas/pkg/hst_calib/stis/defringe.cl>`_.

    Parameters
    ----------
    inspec: str
        Name of input science spectrum
    inflat: str
        Name of input [final] fringe flat, typically derived from `mkfringeflat`
    outspec: str
        Name of output science spectrum
    do_cal: bool
        Perform calstis processing before defringing?  (Default=False)
    darkfile: str or None
        Name of superdark image.  If None, use DARKFILE in main header of input spectrum.
    pixelflat: str or None
        Name of pixel-to-pixel flat.  If None, use PFLTFILE in main header of the inflat.
    do_defringe: bool
        Correct spectrum for fringing?  (Default=True)

    Returns
    -------
    outname: str
        Fully qualified name of output science spectrum

    """
    # These notes are based on the old STSDAS algorithm.

    raise NotImplementedError()


def call_defringe():
    ''' Command line entry point for defringe().
    '''
    import argparse

    parser = argparse.ArgumentParser(
        description='Correct STIS CCD G750L or G750M spectrum for fringing')
    parser.add_argument('inspec', type=str,
        help='Name of input science spectrum, typically derived from mkfringeflat')
    parser.add_argument('inflat', type=str, help='Name of input [final] fringe flat')
    parser.add_argument('outspec', type=str, help='Name of output science spectrum')
    parser.add_argument('--do_cal', '-c', action='store_true',
        help='Perform calstis processing before defringing?')
    parser.add_argument('--darkfile', type=str, default=None,
        help='Name of superdark image. If omitted, use DARKFILE in main header of the '
             'inspec.')
    parser.add_argument('--pixelflat', type=str, default=None,
        help='Name of pixel-to-pixel flat. If omitted, use PFLTFILE in main header of '
             'the inflat.')
    parser.add_argument('--skip_defringe', dest='do_defringe', action='store_false',
        help='Skip correcting the spectrum for fringing?')

    args = vars(parser.parse_args())
    defringe(**args)


if __name__ == '__main__':
    call_defringe()
