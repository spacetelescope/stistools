#! /usr/bin/env python

import os
import numpy as np
import warnings
from astropy.io import fits
from scipy.interpolate import interp1d
from astropy.nddata import utils

from ..r_util import expandFileName
from ..calstis import calstis


# Keyword choices for calstis reduction:
PERFORM = {
    'ALL':   ['DQICORR', 'BLEVCORR', 'BIASCORR', 'DARKCORR', 'FLATCORR', 'CRCORR'],
    'G750M': ['WAVECORR', 'HELCORR', 'X2DCORR'],
    'G750L': [],}
OMIT_G750L = PERFORM['G750M'][:]


def normspflat(inflat, outflat='.', do_cal=True, biasfile=None, darkfile=None,
               pixelflat=None, wavecal=None):
    """Normalize STIS CCD fringe flat.

    Based on PyRAF `stsdas.hst_calib.stis.normspflat task 
    <https://github.com/spacetelescope/stsdas/blob/master/stsdas/pkg/hst_calib/stis/normspflat.cl>`_.

    Parameters
    ----------
    inflat: str
        Name of input fringe flat
    outflat: str
        Name of normalized fringe flat output or directory location.  (Default=".")
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

    # Determine properties of inflat (assumed to be {RAW, CRJ, SX2, X2D, clff}):
    inflat = os.path.abspath(inflat)
    if not os.access(inflat, os.R_OK|os.W_OK):
        raise FileNotFoundError('Cannot access inflat:  {}'.format(inflat))
    hdr0 = fits.getheader(inflat, ext=0)
    rootname = hdr0['ROOTNAME']  # *** Or use something derived from os.path.basename(inflat)? ***
    opt_elem = hdr0['OPT_ELEM'].upper()
    if opt_elem not in ['G750L', 'G750M']:
        raise ValueError('Unsupported opt_elem="{}"!'.format(opt_elem))

    # *** TO DO:  Confirm inflat has >= 2 SCI exts to allow CR-rejection *** line 185 in cl
    # *** TO DO:  Handle NRPTEXP keyword *** line 194 in cl

    # Identify wavecal:
    if opt_elem == 'G750M':
        wavecal = wavecal or expandFileName(hdr0['WAVECAL'])
        if not os.access(wavecal, os.R_OK):
            raise FileNotFoundError('Cannot find G750M wavecal:  {}'.format(wavecal))
    if wavecal in [None, 'N/A', 'n/a', '']:
        wavecal = ''  # Format for non-specified wavecals going into calstis API
    else:
        wavecal = os.path.abspath(wavecal)

    # Handle outflat filename formatting:
    outflat = os.path.abspath(outflat)
    if os.path.isdir(outflat):
        outflat = os.path.join(outflat, '{}_nsp.fits'.format(rootname))
    if not os.access(os.path.dirname(outflat), os.W_OK):
        raise IOError('Do not have permission to write to:  {}'.format(
                      os.path.dirname(outflat)))

    perform = PERFORM['ALL'] + PERFORM[opt_elem]
    if do_cal:
        # Resolve reference files from supplied variables or header values:
        ref_types = {
            'BIASFILE': os.path.abspath(biasfile or expandFileName(hdr0['BIASFILE'])),
            'DARKFILE': os.path.abspath(darkfile or expandFileName(hdr0['DARKFILE'])),
            'PFLTFILE': os.path.abspath(pixelflat or expandFileName(hdr0['PFLTFILE'])),}

        # Populate/repopulate the inflat header accordingly:
        for i, (ref_type, ref) in enumerate(ref_types.items()):
            if not os.access(ref, os.F_OK):
                raise FileNotFoundError('Cannot access reference file:  {}'.format(ref))

            # Handle reference file paths via environment variables:
            ref_var = 'reff{:.0f}'.format(i+1)
            os.environ[ref_var] = os.path.abspath(os.path.dirname(ref)) + os.sep
            # Keep $oref where it's the same:
            if os.path.normpath(os.environ[ref_var]) == os.path.normpath(os.environ['oref']):
                ref_var = 'oref'
                #os.environ['oref'] = os.path.abspath(os.environ['oref'])  # for when we cd
            with fits.open(inflat, 'update') as f:
                if ref_var == 'oref':
                    f[0].header[ref_type] = '{}${}'.format(ref_var, os.path.basename(ref))
                else:
                    f[0].header[ref_type] = '${}/{}'.format(ref_var, os.path.basename(ref))

        # Update calibration flags prior to calling calstis:
        with fits.open(inflat, 'update') as f:
            for keyword in perform:
                if f[0].header[keyword] != 'COMPLETE':
                    f[0].header[keyword] = 'PERFORM'
            if opt_elem == 'G750M':
                output_suffix = 'sx2'  # geometrically-corrected
            elif opt_elem == 'G750L':
                output_suffix = 'crj'  # not geometrically-corrected
                for keyword in OMIT_G750L:
                    if f[0].header[keyword] == 'PERFORM':
                        f[0].header[keyword] = 'OMIT'

        # Call calstis:
        old_cwd = os.getcwd()  # outflat and ref_vars are abs paths
        try:
            os.chdir(os.path.dirname(inflat))  # Must be in same directory to pick up EPC file
            trailer = os.path.join(os.path.dirname(outflat),
                                   '{}_calstis.log'.format(rootname))
            outroot = os.path.dirname(outflat) + os.sep
            outname = os.path.join(os.path.dirname(outflat),
                                   '{}_{}.fits'.format(rootname, output_suffix))
            # Remove files from previous iterations of this script:
            for old_file in [outname, trailer]:
                if os.access(old_file, os.F_OK):
                    os.remove(old_file)
            res = calstis(os.path.basename(inflat), wavecal=wavecal, outroot=outroot,
                          trailer=trailer)
            if res != 0:
                raise Exception('CalSTIS returned non-zero code:  {}\n' \
                                'See log for more details:  {}'.format(res, trailer))
            print ('File written:  {}'.format(outname))
        finally:
            os.chdir(old_cwd)
    else:  # not do_cal
        outname = inflat

        # Check if user-supplied inflat has all the recommended calibrations performed:
        keyword_warnings = []
        for keyword in perform:
            if hdr0[keyword] != 'COMPLETE':
                keyword_warnings.append(keyword)
        if keyword_warnings:
            warnings.warn('These calibration steps should be COMPLETE:\n{}'.format(
                ', '.join(keyword_warnings)))

    # Read in the calibrated flat data:
    data = fits.getdata(outname, ext=1)

    # Get shape information, CL gets this from the header, but simpler to get it from the data?
    numrows, numcols = np.shape(data)

    # Generate clff file, this may be unneccesary in python
    if do_cal:
        flat_used = rootname
    else:
        flat_used = inflat
        if opt_elem != "G750M":
            print("It is assumed that pixel-to-pixel flat fielding has already been performed.")

    if opt_elem == "G750M":
        with fits.open(flat_used+"_sx2.fits") as hdulist:
            # generate a _tmp.fits file off the _sx2 file divided by itself, this is just a ones array

            hdulist[1].data = np.ones(np.shape(hdulist[1].data))
            hdulist.writeto(str(outflat.split("_")[0]) + "_tmp.fits") # Do we really need to generate intermediate products in python?
    else:
        with fits.open(flat_used + "_crj.fits") as hdulist:
            # generate a clff file, which is just a copy of the crj file?
            hdulist.writeto(rootname + "_clff.fits")
            # generate a _tmp.fits file off the crj file divided by itself, this is just a ones array
            hdulist[1].data = np.ones(np.shape(hdulist[1].data))
            hdulist.writeto(str(outflat.split("_")[0]) + "_tmp.fits")  # Do we really need to generate intermediate products in python?


    # if do_cal
    # # if G750M generate an out _tmp.fits file off the _sx2 file divided by itself (?)
    # # else if theres a clff file delete it and make a new clff file thats the crj divided by the pfltfile
    # else (if not do_cal)
    # # if G750M make a _tmp file of the inflat divided by itself, this might just be a 1s array with all the hdr info
    # # else (if not G750M) assume that pixel flat fielding has been done and make a tmp file that is the inflat divided
    # # by itself and make a clff file



    # Do a line-by-line cubic spline fit to the fringe flat to remove the lamp function...

    #blkavg can be done using astropy.nddata.utils.block_reduce

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
    parser.add_argument('--outflat', '-o', type=str, default='.',
        help='Name of normalized fringe flat output or directory location (default=".")')
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
