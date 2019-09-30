#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits

from ..r_util import expandFileName


def mkfringeflat(inspec, inflat, outflat, do_shift=True,
                 shift_range=[-0.5, 0.5], shift_step=0.1,
                 do_scale=True, scale_range=[0.8, 1.2], scale_step=0.04,
                 extrloc=None, extrsize=None, opti_spreg=None, rmsregion=None):
    """Match fringes in STIS fringe flat to those in science data.

    Based on PyRAF `stsdas.hst_calib.stis.mkfringeflat task
    <https://github.com/spacetelescope/stsdas/blob/master/stsdas/pkg/hst_calib/stis/mkfringeflat.cl>`_.

    Parameters
    ----------
    inspec: str
        Name of input science spectrum
    inflat: str
        Name of input normalized fringe flat from `normspflat`
    outflat: str
        Name of output [final] fringe flat
    do_shift: bool
        Shift flat to match fringes in spectrum?  (Default=True)
    shift_range: [float, float]
        Range for shift determination.  (Default=[-0.5, 0.5])
    shift_step: float or None
        Step size used when calculating shifts.  (Default=0.1)
    do_scale: bool
        Scale fringe amplitude to match science spectrum?  (Default=True)
    scale_range: [float, float]
        Range for scale determination.  (Default=[0.8, 1.2])
    scale_step: float
        Step size used when calculating scale factors.  (Default=0.04)
    extrloc: float or None
        Central line (row) to be extracted from 2-D spectrum (zero-indexed).  If None, ...
    extrsize: float or None
        Number of lines to be extracted from 2-D spectrum.  If None, ...
    opti_spreg: [int, int] or None
        Spectral range used in normalizing the spectrum, as specified via zero-indexed
        pixel numbers.  If None, ...
    rmsregion: [int, int] or None
        Spectral range used in measuring RMSs, as specified via zero-indexed pixel
        numbers.  If None, ...

    Returns
    -------
    outname: str
        Fully qualified name of the outflat

    """
    # These notes are based on the old STSDAS algorithm.

    # check input flags:  at least one of do_shift & do_scale should be set
    # check whether inflat has data in ext=0 or ext=[SCI,1] --> op1file
    # if not specified, determine EXTRSIZE using APERTURE of fringe flat
    # if EXTRLOC was specified, validate it
    # if not specified, determine EXTRLOC:
    #    use row with max counts in the spectrum from the middle 60% of cols
    #    refine peak with parabola fit near max row
    # determine actual rows from EXTRSIZE & EXTRLOC
    # interpret opti_spreg region wrt binning
    # if not specified, calculate default rmsregion

    # if do_shift:
    #    determine discreet shifts from shift_range & shift_step
    #    shift inflat fractional-pixel amount via linear interpolation in x
    #    calculate inspec / shifted_inflat
    #    if G750M:
    #        noao.twodspec.longslit.response with cubic spline with 1 piece
    #    elif G750L:
    #        noao.twodspec.longslit.response with cubic spline with 15 pieces
    #    collapse response down to 1d using extraction region
    #    calculate stddev/mean within rmsregion
    #    find shift that minimizes stddev/mean:
    #        if min shift is near edge of range (within 2 steps):
    #            issue warning
    #            return edge shift value
    #        else:
    #            return weighted shift value from 5 locations nearest min
    #    apply the best shift --> inflat_sh
    #    mark output file with keyword SHIFTED

    # if do_scale:
    #    parse opti_spreg for range (not necessary with list input here)
    #    modify NEXTEND keyword (?)
    #    determine range of scale factors to try
    #    scale inflat (shifted?)
    #    if G750M:
    #        noao.twodspec.longslit.response with cubic spline with 1 piece
    #    elif G750L:
    #        noao.twodspec.longslit.response with cubic spline with 15 pieces
    #    collapse response down to 1d using extraction region
    #    calculate stddev/mean within rmsregion
    #    find scale that minimizes stddev/mean:
    #        if min scale is near edge of range (within 2 steps):
    #            issue warning
    #            return edge scale value
    #        else:
    #            return weighted scale value from 5 locations nearest min
    #    apply the best scale --> theoutflat
    #    mark output file with keyword SCALED (not actually done)

    # save outflat

    raise NotImplementedError()


def call_mkfringeflat():
    """Command line entry point for mkfringeflat().
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Match fringes in STIS fringe flat to those in science data')
    parser.add_argument('inspec', type=str, help='Name of input science spectrum')
    parser.add_argument('inflat', type=str,
        help='Name of input normalized fringe flat from normspflat')
    parser.add_argument('outflat', type=str, help='Name of output [final] fringe flat')
    parser.add_argument('--skip_shift', dest='do_shift', action='store_false',
        help='Skip shifting the flat to match fringes in spectrum?')
    parser.add_argument('--shift_range', type=float, default=[-0.5, 0.5], nargs=2,
        metavar='FLOAT', help='Range for shift determination (default=[-0.5, 0.5])')
    parser.add_argument('--shift_step', type=float, default=0.1, metavar='FLOAT',
        help='Step size used when calculating shifts (default=0.1)')
    parser.add_argument('--skip_scale', dest='do_scale', action='store_false',
        help='Skip scaling the fringe amplitude to match science spectrum?')
    parser.add_argument('--scale_range', type=float, default=[0.8, 1.2], nargs=2,
        metavar='FLOAT', help='Range for scale determination (default=[0.8, 1.2])')
    parser.add_argument('--scale_step', type=float, default=0.04, metavar='FLOAT',
        help='Step size used when calculating scale factors (default=0.04)')
    parser.add_argument('--extrloc', '-l', type=float, default=None, metavar='FLOAT',
        help='Central line (row) to be extracted from 2-D spectrum (zero-indexed). '
             'If omitted, ...')
    parser.add_argument('--extrsize', '-s', type=float, default=None, metavar='INT',
        help='Number of lines to be extracted from 2-D spectrum. If omitted, ...')
    parser.add_argument('--opti_spreg', type=int, nargs=2, default=None, metavar='INT',
        help='Spectral range used in normalizing the spectrum, as specified via '
             'zero-indexed pixel numbers. If omitted, ...')
    parser.add_argument('--rmsregion', type=int, nargs=2, default=None, metavar='INT',
        help='Spectral range used in measuring RMSs, as specified via zero-indexed pixel '
             'numbers. If omitted, ...')

    args = vars(parser.parse_args())
    mkfringeflat(**args)


if __name__ == '__main__':
    call_mkfringeflat()
