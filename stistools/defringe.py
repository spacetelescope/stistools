#! /usr/bin/env python

import argparse
from astropy.io import fits
import datetime
import os
import textwrap
import re

import numpy as np

from stistools.r_util import expandFileName

#   4  bad detector pixel or beyond aperture
#   8  data masked by occulting bar
# 512  bad pixel in reference file
sdqflags = 4 + 8 + 512


def defringe(science_file, fringe_flat, output=None,
             overwrite=True, verbose=True):
    # dark_file=None, pixel_flat=None
    """Defringe by dividing the science spectrum by the fringe flat

    Parameters
    ----------
    science_file : str
        The name (including directory) of the input science file.

    fringe_flat : str
        The name (including directory) of the input fringe flat file.  This
        is the output from `mkfringeflat`.

    output : str or None
        If `output` is None, the name of the output file will be constructed
        by appending "_defringe" to the root of the input file name.
        If `output` is not None, this is the name (not including directory)
        of the output file.  The output will be written to the same directory
        as the input science file.

    overwrite : bool
        If there is an existing file with the same name as the output name,
        the existing file will be overwritten if `overwrite` is True (the
        default is True).
        If the input and output names are the same, a RuntimeError will be
        raised, rather than modifying the input in-place.

    verbose : bool
        If True (the default), print more info.

    Returns
    -------
    output_filename : str
        The name of the output file, including directory.
    """

    if science_file.endswith("_raw.fits"):
        print('Warning:  The science file name ends with "_raw.fits".  If this '
              'really is raw data, this script will fail.', flush=True)

    # Define new filename:
    science_file = os.path.normpath(expandFileName(science_file))  # Expand IRAF and UNIX $VARS
    sci_dir, sci_filename = os.path.split(science_file)
    if output is None:
        sci_root = re.split(r'\.fits.*', sci_filename,
                            flags=re.IGNORECASE)[0].rsplit('_',1)[0]
        output_filename = os.path.join(sci_dir, sci_root + '_defringe.fits')
    else:
        output_filename = os.path.join(sci_dir, output)
    if science_file == output_filename:
        raise RuntimeError('The input and output file names cannot be the same.')

    # Get the data from the fringe flat file:
    fringe_flat = os.path.normpath(expandFileName(fringe_flat))  # Expand IRAF and UNIX $VARS
    with fits.open(fringe_flat) as fringe_hdu:
        # We assume that either the fringe flat data are in the primary HDU,
        # or there is one imset containing the data.
        if len(fringe_hdu) == 1:
            fringe_data = fringe_hdu[0].data
            fringe_dq   = None
            if verbose:
                print('Fringe flat data were read from the primary HDU')
        else:
            fringe_data = fringe_hdu[("sci",1)].data
            fringe_dq   = fringe_hdu[("dq",1)].data
            if verbose:
                if fringe_dq is None:
                    print('Fringe flat data were read from the first imset')
                else:
                    print('Fringe flat data and DQ were read from the first imset')
    if fringe_data is None:
        raise RuntimeError('There is no data in the fringe flat.')

    # Since we're going to divide by fringe_data, make sure there aren't any
    # pixels where it's zero.  There shouldn't be any negative values, either.
    fringe_mask = (fringe_data <= 0.)
    n_not_positive = fringe_mask.sum()
    if n_not_positive > 0 and verbose:
        print('{} pixels in the fringe flat were less than or equal to 0'
              .format(n_not_positive))
    if fringe_dq is not None:
        temp = np.bitwise_and(fringe_dq, sdqflags)
        fringe_dq_mask = (temp > 0)
        del temp
        n_fringe_dq = fringe_dq_mask.sum()
        if n_fringe_dq > 0:
            fringe_mask = np.logical_or(fringe_mask, fringe_dq_mask)
            if verbose:
                print('{} pixels were flagged as bad in the fringe flat DQ array'
                      .format(n_fringe_dq))
    # For pixels that are bad in the fringe flat, we will not make any change
    # to the science data.
    # Note:  Save fringe_mask and n_fringe_mask for later.
    n_fringe_mask = fringe_mask.sum()
    if n_fringe_mask > 0:
        fringe_data[fringe_mask] = 1.

    # Correct the data in the science file:
    with fits.open(science_file) as science_hdu:
        n_hdu = len(science_hdu)
        # Get a list of all EXTVER values in the science file
        extver = []
        for hdunum in range(1, n_hdu):
            hdr = science_hdu[hdunum].header
            if 'extver' in hdr:
                extver.append(int(hdr['extver']))
        # Remove duplicates, and sort
        imsets = set(extver)

        # Apply the fringe flat to all image sets.
        for extver in imsets:
            try:
                science_data = science_hdu[('sci', extver)].data
            except KeyError:
                print('Warning:  HDU ("SCI", {}) not found'.format(extver))
                science_data = None
            try:
                science_dq = science_hdu[('dq', extver)].data
            except KeyError:
                print('Warning:  HDU ("DQ", {}) not found'.format(extver))
                science_dq = None
            try:
                # Note that science_err will be None if this HDU has no data array.
                science_err = science_hdu[('err', extver)].data
            except KeyError:
                print('Warning:  HDU ("ERR", {}) not found'.format(extver))
                science_err = None

            # Divide science data and error arrays by fringe flat:
            if science_data is not None:
                science_data /= fringe_data
            if science_err is not None:
                science_err  /= fringe_data

            # Update the DQ array in the science file.
            if science_dq is None:
                science_dq = np.zeros(science_data.shape, dtype=np.int16)
            if n_fringe_mask > 0:
                # Flag pixels that had fringe flat data <= 0 or were flagged
                # as bad in the fringe flat DQ array.
                science_dq[fringe_mask] |= 512          # bad pixel in ref file
            if fringe_dq is not None:
                # Combine fringe flat DQ with science DQ.
                science_dq = np.bitwise_or(fringe_dq, science_dq)
            science_hdu[('dq', extver)].data = science_dq

            if verbose:
                print('Imset {} done'.format(extver))

        # Update primary header:
        science_hdu[0].header.add_history(' ')
        science_hdu[0].header.add_history('DEFRINGING complete ...')
        science_hdu[0].header.add_history('  reference fringe flat {}'.format(os.path.basename(fringe_flat)))
        science_hdu[0].header.add_history('  {}'.format(str(datetime.datetime.now()).rsplit('.')[0]))
        science_hdu[0].header.add_history('  science and error arrays divided by fringe flat.')

        # Write to a new FITS file
        # Remove old version, if it exists:
        if os.path.exists(output_filename) and overwrite:
            print('Removing and recreating {}'.format(output_filename))
            os.remove(output_filename)
        science_hdu.writeto(output_filename)
        print('Defringed science data saved to {}'.format(output_filename))

    return output_filename


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
        Steps
        --------------------------------------------------------------
        1. Get the data from the science and fringe flat files
        2. Replace zero or negative values in the fringe flat with 1
        3. Divide science data and error arrays by fringe flat
        4. Flag DQ with 512 where the fringe flat was <= 0
        5. Added HISTORY entries in the header
        6. Write to a new fits file

        Output
        ----------------------------------------------------------------
        Outputs a file that has the same rootname as the input science file,
        but with '_defringe.fits' at the end. This is the final defringed data.
        '''), description='Script to calibrate science data in preparation for defringing')

    parser.add_argument('science',
                        type=str,
                        help='file containing the calibrated science data. This should be '
                             'the output from prepspec')
    parser.add_argument('fringeflat',
                        type=str,
                        help='file containing the fringe flat. This should be the output '
                             'from mkfringeflat')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    defringe(args.science, args.fringeflat)
