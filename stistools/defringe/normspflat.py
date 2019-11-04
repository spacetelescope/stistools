#! /usr/bin/env python

import os
import numpy as np
import warnings
from astropy.io import fits
from astropy.nddata import utils

from ..r_util import expandFileName
from ..calstis import calstis
from ..fit1d import fit1d


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
    numrows, numcols = np.shape(data)

    # Do a line-by-line cubic spline fit to the fringe flat to remove the lamp function...

    with fits.open(outname) as hdulist:

        aperture = hdulist[0].header['APERTURE']
        bincols = hdulist[0].header['BINAXIS1']
        binrows = hdulist[0].header['BINAXIS2']
        cenwave = hdulist[0].header['CENWAVE']

        # If short-slit aperture (does not start with "52X"):
        if aperture[0:3] != "52X":
            flatdata = hdulist[1].data
            #    Find the row with max counts in the short-slit fringe flat
            #    Search between the middle 10% of rows and the middle 60% of columns
            t_row = int(0.45 * numrows)
            b_row = int(0.55 * numrows)
            #     Uses central 60% of column
            l_col = int(0.2 * numcols)
            r_col = int(0.8 * numcols)

            row_avgs = np.array([np.average(row) for row in flatdata[t_row:b_row, l_col:r_col]])
            max_row_idx = np.where(row_avgs == np.max(abs(row_avgs)))[0][0] + t_row  # CL does an absolute value here
            max_row = flatdata[max_row_idx]
            #    Does some flux-filtering -- I think this is for determining the max rows


        # Set rows (startrow, lastrow) to be fit according to aperture name (and possibly OPT_ELEM):
        # '0.3X0.09', '0.2X0.06', '52X...'

        if aperture == "0.3X0.09":
            startrow = max_row_idx - int(4./binrows+0.25)
            lastrow = max_row_idx + int(4./binrows+0.25)
        elif aperture == "0.2X0.06":
            startrow = max_row_idx - int(3. / binrows + 0.25)
            lastrow = max_row_idx + int(3. / binrows + 0.25)
        elif aperture[0:3] != "52X":
            print("ERROR: not able to understand APERTURE keyword")
            return
        elif aperture == "G750M":
            startrow = int(92./binrows) + 1
            lastrow = startrow + int(1024./binrows) - 1
        else:
            startrow = 1
            lastrow = numrows
        print(startrow, lastrow, max_row_idx)

        # Details of spline fit determined according to OPT_ELEM + CENWAVE.
        # G750M (various CENWAVEs), G750L (i.e. CENWAVE == 7751; various binning), other (never used?)
        if opt_elem == "G750M":
            fitted = []  # This probably needs to be padded with ones to match the input shape

            if cenwave < 9800:
                for row in data:
                    l_row_data = row[0:86]
                    r_row_data = row[1109:]
                    fit_data = row[86:1109]
                    xrange = np.arange(0, len(fit_data), 1.)
                    spl = fit1d(xrange, fit_data, naverage=2, function="spline3",
                                order=1, low_reject=5.0, high_reject=5.0, niterate=2)
                    row_fit = spl(xrange)
                    fitted_row = np.array(
                        list(np.ones(len(l_row_data))) + list(row_fit) + list(np.ones(len(r_row_data))))
                    fitted.append(fitted_row)
                pass

            elif cenwave == 9851:
                for row in data:
                    l_row_data = row[0:86]
                    r_row_data = row[1109:]
                    fit_data = row[86:1109]
                    xrange = np.arange(0, len(fit_data), 1.)
                    spl = fit1d(xrange, fit_data, naverage=2, function="spline1",
                                order=2, low_reject=5.0, high_reject=5.0, niterate=2)
                    row_fit = spl(xrange)
                    fitted_row = np.array(
                        list(np.ones(len(l_row_data))) + list(row_fit) + list(np.ones(len(r_row_data))))
                    fitted.append(fitted_row)
            else:
                for row in data:
                    l_row_data = row[0:86]
                    r_row_data = row[1109:]
                    fit_data = row[86:1109]
                    xrange = np.arange(0, len(fit_data), 1.)
                    spl = fit1d(xrange, fit_data, naverage=2, function="spline3",
                                order=2, low_reject=5.0, high_reject=5.0, niterate=2)
                    row_fit = spl(xrange)
                    fitted_row = np.array(
                        list(np.ones(len(l_row_data))) + list(row_fit) + list(np.ones(len(r_row_data))))
                    fitted.append(fitted_row)

            # Write to the output file
            hdulist[1].data = np.array(fitted)
            hdulist.writeto(str(outflat.split(".")[0]) + ".fits")
            return

        elif cenwave == 7751:  # G750L

            if bincols == 1:
                startcol = 591
                endcol = 640
                highorder = 60
                loworder = 12
            elif bincols == 2:
                startcol = 295
                endcol = 320
                highorder = 50
                loworder = 12
            elif bincols == 4:
                startcol = 145
                endcol = 160
                highorder = 50
                loworder = 12

            # Fit both the high order and low order splines
            fitted_highorder = []
            fitted_loworder = []
            for row in data:
                fit_data = row[:]
                xrange = np.arange(0, len(fit_data), 1.)

                # High Order Fit
                spl = fit1d(xrange, fit_data, naverage=2, function="spline3",
                            order=highorder, low_reject=5.0, high_reject=5.0, niterate=2)
                fitted_row = spl(xrange)
                fitted_highorder.append(fitted_row)

                # Low Order Fit
                spl = fit1d(xrange, fit_data, naverage=2, function="spline3",
                            order=loworder, low_reject=5.0, high_reject=5.0, niterate=2)
                fitted_row = spl(xrange)
                fitted_loworder.append(fitted_row)

            fitted_highorder = np.array(fitted_highorder)
            print(fitted_highorder[512])
            fitted_loworder = np.array(fitted_loworder)
            # Divide both spline fits off the science data
            resp_highorder = data/fitted_highorder
            print(resp_highorder[512])
            resp_loworder = data/fitted_loworder

            # Get absolute difference and ratio of the two response arrays
            resp_absdiff = abs(resp_highorder - resp_loworder)
            resp_ratio = resp_highorder/resp_loworder

            # Iterate through the rows from startrow to lastrow
            fitted = []
            for row_idx in np.arange(0, startrow, 1):
                fitted.append(list(np.ones(numcols)))

            for row_idx in np.arange(startrow, lastrow, 1):
                # Find the min pixel location between startcol and endcol
                min_idx = np.argmin(resp_absdiff[row_idx][startcol:endcol]) + startcol  # the idx in the full array

                # Use the ratio between the two fits at the minimum point as the scale factor between the two fits
                scale_factor = resp_ratio[row_idx][min_idx]
                print(scale_factor)
                # Generate the flat using the high order flat for any column left of the min pixel and the low order
                # flat times the scale factor for anything right of the min pixel
                highorder_segment = resp_highorder[row_idx][:min_idx+1]
                loworder_segment = resp_loworder[row_idx][min_idx+1:] * scale_factor

                stitched_row = np.array(list(highorder_segment) + list(loworder_segment))
                fitted.append(stitched_row)

            for row_idx in np.arange(lastrow, numcols, 1):
                fitted.append(list(np.ones(numrows)))

            fitted = np.array(fitted)

            # Write to the output file
            hdulist[1].data = np.array(fitted)
            hdulist.writeto(str(outflat.split(".")[0]) + ".fits")
            return

        else:
            fitted = []
            for row in data:
                xrange = np.arange(0, len(row), 1.)
                spl = fit1d(xrange, row, naverage=2, function="spline3",
                            order=20, low_reject=3.0, high_reject=3.0, niterate=2)
                row_fit = spl(xrange)
                fitted.append(row_fit)

            # Write to the output file
            hdulist[1].data = np.array(fitted)
            hdulist.writeto(str(outflat.split(".")[0]) + ".fits")
            return



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
