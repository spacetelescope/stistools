#! /usr/bin/env python

from astropy.io import fits
from astropy.nddata.blocks import block_reduce
import numpy as np
import math
import os
from scipy.ndimage import shift
from ._findloc import find_loc
from ._response import response

__version__ = 0.1

def mkfringeflat(inspec, inflat, outflat, do_shift=True, beg_shift=-0.5, end_shift=0.5,
                 shift_step=0.1, do_scale=True, beg_scale=0.8, end_scale=1.2, scale_step=0.04,
                 extrloc=None, extrsize=None, opti_spreg=None, rms_region=None):

    """Takes an input science spectrum and a fringe flat that has been normalized using
    the task `normspflat`.  The fringe flat is shifted and scaled to produce the minimum
    RMS when divided into the science data.

    Based on the PyRAF `stsdas.hst_calib.stis.mkfringeflat` task.

    In `mkfringeflat`, the user can specify a range of shifts and scales for the routine
    to test creating an optimal fringe flat. `mkfringeflat` will go through the shift and
    scale dimensions separately and calculate the RMS using the following steps:

    1. For each shift step, apply the shift to the input flat field

    2. Divide the science data by the shifted flat

    3. Divide out the large-scale SED from the science image using a spline fit in order
       to isolate the fringing pattern (this is called the response image)

    4. Sum the response image along the columns within the RMS region

    5. Calculate the mean and standard deviation of the summed columns of the response
       image

    6. The RMS value for that shift is given by the standard deviation divided by the mean
       found in step 5

    7. Fit the RMS values with a quadratic polynomial weighted by the inverse RMS to find
       the optimal RMS value

    8. Apply the best shift determined in step 7 to the data and repeat steps 1-7 with the
       scale values to find the best scaling


    The RMS values are printed out for each scale and shift but the final best shift and
    best scale values do not necessarily correspond to the printed values. This is
    because the routine is calculating the RMS values based on a fit of the data at each
    scale and shift, rather than being calculated at each discrete step.

    Parameters
    ----------

    inspec: str
        Name of input science spectrum datafile

    inflat: str
        Name of input fringe flat file (usually the output from `normspflat`)

    outflat: str
        Name of output fringe flat to be used in the defringe task

    do_shift: bool
        Controls whether the shift between fringe flat and science data is
        to be calculated

    beg_shift: float
        Initial shift to apply to fringe flat

    end_shift: float
        Final shift to apply to fringe flat

    shift_step: float
        Step-size between shifts to be applied to fringe flat

    do_scale: bool
        Controls whether the scaling between fringe flat and science
        data is to be calculated

    beg_scale: float
        Initial scaling to apply to fringe flat

    end_scale: float
        Final scaling to appply to fringe flat

    scale_step: float
        Step-size between scaling values to be applied to fringe flat

    extrloc: float or None
        Extraction location.  If set to None, this will be calculated by
        parabolic interpolation of the peak of the cross-dispersion
        spectral sum

    extrsize: float or None
        Extraction size in pixels.  If set to None, this will be set to a
        reasonable value by this routine

    opti_spreg: list or array-like or None
        A list or array representing the section to be used in normalizing the spectrum
        of the science target before it is divided by the shifted/scaled fringe flat.
        If set to None, a reasonable range is chosen by this routine.  Should be
        specified like a Python slice, zero indexed.

    rms_region: list or array-like or None
        A list or array representing the section to be used in the rms calculation.  If set
        to None, a reasonable range is chosen by this routine.  Should be specified
        like a Python slice, zero indexed.
    """

    print("mkfringeflat.py version {}".format(__version__))
    print(" - matching fringes in a flatfield to those in science data")

    sci_hdulist = fits.open(inspec)
    flt_hdulist = fits.open(inflat)

    scidata = sci_hdulist[1].data
    sci_prihdr = sci_hdulist[0].header
    sci_hdr = sci_hdulist[1].header

    if len(flt_hdulist) == 1:
        fltdata = flt_hdulist[0].data
        flt_hdr = flt_hdulist[0].header
    else:
        fltdata = flt_hdulist[1].data
        flt_hdr = flt_hdulist[1].header

    flt_prihdr = flt_hdulist[0].header

    try:
        shifted = flt_prihdr['shifted']
    except KeyError:
        shifted = ""
    if "YE" in shifted:
        print(" ")
        print(" NOTE: Input flat was already shifted")

    opt_elem = sci_prihdr['opt_elem']

    nrows, ncols = scidata.shape
    bincols = sci_prihdr['binaxis1']
    binlines = sci_prihdr['binaxis2']

    fltbincols = flt_prihdr['binaxis1']
    fltbinlines = flt_prihdr['binaxis2']

    ltv1 = sci_hdr['ltv1']

    if opt_elem == "G750M":

        centera2 = sci_prihdr['centera2']
        sizaxis2 = sci_prihdr['sizaxis2']

        ltv2 = 1 - centera2 + sizaxis2/2.0

        flt_centera2 = flt_prihdr['centera2']
        flt_sizaxis2 = flt_prihdr['sizaxis2']

        flt_ltv2 = 1 - flt_centera2 + flt_sizaxis2/2.0

    else:

        ltv2 = sci_hdr['ltv2']
        flt_ltv2 = flt_hdr['ltv2']

    sax0 = round(-ltv1) + 1
    sax1 = sax0 + ncols - 1

    say0 = int(round(flt_ltv2 - ltv2)) + 1
    say1 = say0 + nrows - 1

    #
    # convert IRAF section limits to Python slice indices

    shiftrowstart = say0 - 1
    shiftrowstop = say1
    shiftcolstart = sax0 - 1
    shiftcolstop = sax1

    aperture = sci_prihdr['aperture']

    if extrsize is None:
        if "0.3X0.09" in aperture:
            numextrows = int(round(9.0 / binlines + 0.25))
        elif "0.2X0.06" in aperture:
            numextrows = int(round(7.0 / binlines + 0.25))
        elif "52X" not in aperture:
            print("ERROR: not able to understand APERTURE keyword")
            return
        else:
            numextrows = 11/binlines
        apername = " [Aperture: " + aperture + "]"
    else:
        numextrows = extrsize
        apername = ""

    if extrloc is None:
        maxrow = int(round(find_loc(scidata)))
    else:
        maxrow = extrloc

    print(" Extraction center: row {}".format(maxrow))
    print("   Extraction size: {} pixels {}".format(numextrows, apername))
    fline = maxrow - int(round((numextrows - 0.49999)/2.0))
    lline = maxrow + int(round((numextrows - 0.49999)/2.0)) + 1

    if opti_spreg is None:
        if "G750M" in opt_elem:
            colstart = int(83/bincols)
            colstop = int(1106/bincols)
        else:
            colstart = int(5/bincols) - 1
            colstop = int(1020/bincols)
    else:
        colstart, colstop = opti_spreg
        colstart, colstop = int(colstart), int(colstop)

    if rms_region is None:
        rms_start = int(725/bincols) - 1
        rms_stop = int(900/bincols)
    else:
        rms_start, rms_stop = rms_region
        rms_start, rms_stop = int(rms_start), int(rms_stop)

    print("Range to be normalized: [{}:{},{}:{}]".format(fline, lline, colstart, colstop))

    shifted_flat = None

    if do_shift:
        print("")
        print("Determining best shift for fringe flat")
        print("")

        expo = int(round(math.log10(shift_step)-0.49999999))

        fshift = int(round(beg_shift/(10**(expo))))
        lshift = int(round(end_shift/(10**(expo))))
        step = int(round(shift_step/(10**(expo))))
        flt_blk = block_reduce(fltdata, (binlines/fltbinlines, bincols/fltbincols),
                                func=np.mean)
        nshifts = (lshift - fshift)//step + 1

        rmsvalues = np.zeros(nshifts)
        current_shift = np.zeros(nshifts)

        for i in range(nshifts):
            current_shift[i] = beg_shift + i*shift_step

            shifted_flat = shift(flt_blk[shiftrowstart:shiftrowstop, shiftcolstart:shiftcolstop],
                                 (0, current_shift[i]), order=1, mode='nearest')

            star_cont_shift = scidata / shifted_flat

            if "G750M" in opt_elem:
                norm_star_cont_shift = response(star_cont_shift[fline:lline,colstart:colstop],
                                                star_cont_shift[fline:lline,colstart:colstop],
                                                threshold=1.0e-20,
                                                function="spline3",
                                                sample="*",
                                                naverage=2,
                                                order=1,
                                                low_reject=3.0,
                                                high_reject=3.0,
                                                niterate=2,
                                                grow=0.0)
            elif "G750L" in opt_elem:
                norm_star_cont_shift = response(star_cont_shift[fline:lline,colstart:colstop],
                                                star_cont_shift[fline:lline,colstart:colstop],
                                                threshold=1.0e-20,
                                                function="spline3",
                                                sample="*",
                                                naverage=2,
                                                order=15,
                                                low_reject=3.0,
                                                high_reject=3.0,
                                                niterate=2,
                                                grow=0.0)
            summed_line = norm_star_cont_shift.sum(axis=0)
            mean = np.mean(summed_line[rms_start:rms_stop], dtype=np.float64)
            sigma = np.std(summed_line[rms_start:rms_stop], dtype=np.float64)
            rmsvalues[i] = sigma/mean
            print("shift = {:10.3f}, rms = {:8.4f}".format(current_shift[i], rmsvalues[i]))

        #
        # Determine shift that delivers the best RMS by an inverse rms weighted average
        weight = 1.0/(rmsvalues*rmsvalues)
        rownum = range(nshifts)
        weighted_rownum = weight*rownum
        minrms = rmsvalues.min()
        min_row = np.where(rmsvalues == minrms)[0][0]

        if min_row == 0 or min_row == nshifts - 1:
            print(" ")
            print("WARNING: Best shift found on the edge of the specified shift range.")
            print("You are advised to try again after adjusting the shift range accordingly")
            theshift = current_shift[min_row]

        else:
            if min_row >= 2 and min_row <= nshifts - 3:
                first_row = min_row - 2
                last_row = min_row + 3
            elif min_row == 1 or min_row == nshifts - 2:
                first_row = min_row - 1
                last_row = min_row + 2
            w_shift_av = weighted_rownum[first_row:last_row].sum()
            weight_av = weight[first_row:last_row].sum()
            theshift = w_shift_av/weight_av
            theshift = theshift*shift_step + beg_shift

        print(" ")
        print(" Best shift : {:10.3f} pixels".format(theshift))

        # Apply the best shift and create output array
        flt_blk = block_reduce(fltdata, (binlines/fltbinlines, bincols/fltbincols),
                                func=np.mean)
        shifted_flat = shift(flt_blk[shiftrowstart:shiftrowstop, shiftcolstart:shiftcolstop],
                             (0, theshift), order=1, mode='nearest')

        # Write out file
        fitspos = inflat.find('.fits')
        output_filename = inflat[:fitspos] + '_sh.fits'
        if os.path.exists(output_filename):
            os.remove(output_filename)
        fits.writeto(output_filename, data=shifted_flat)
        fits.open(output_filename, mode='update')[0].header['shifted'] = 'YES'

        print(" Shifted flat : {}".format(output_filename))
        print("                (Can be used as input flat for next iteration)")

    if do_scale:
        print("")
        print("Determining best scaling of amplitude of fringes in flat")
        print("")


        ilow = int(round(beg_scale*100.0))
        iupp = int(round(end_scale*100.0))
        istep = int(round(scale_step*100.0))

        fltdata = get_flat_data(inflat, shifted_flat)

        flat_mean = fltdata[fline:lline, colstart:colstop].mean()

        nscales = int(round((end_scale - beg_scale)/scale_step)) + 1

        rmsvalues = np.zeros(nscales)
        current_scale = np.zeros(nscales)


        for i in range(nscales):

            current_scale[i] = beg_scale + i*scale_step

            flat_scaled = fltdata.copy()
            flat_scaled[:, colstart+1:colstop-1] = (flat_scaled[:, colstart+1:colstop-1] - flat_mean) * current_scale[i] + \
                                                   flat_mean
            star_cont_scale = scidata / flat_scaled

            if "G750M" in opt_elem:
                norm_star_cont_scale = response(star_cont_scale[fline:lline,colstart:colstop],
                                                star_cont_scale[fline:lline,colstart:colstop],
                                                threshold=1.0e-20,
                                                function="spline3",
                                                sample="*",
                                                naverage=2,
                                                order=1,
                                                low_reject=3.0,
                                                high_reject=3.0,
                                                niterate=2,
                                                grow=0.0)
            elif "G750L" in opt_elem:
                norm_star_cont_scale = response(star_cont_scale[fline:lline,colstart:colstop],
                                                star_cont_scale[fline:lline,colstart:colstop],
                                                threshold=1.0e-20,
                                                function="spline3",
                                                sample="*",
                                                naverage=2,
                                                order=15,
                                                low_reject=3.0,
                                                high_reject=3.0,
                                                niterate=2,
                                                grow=0.0)
            summed_line = norm_star_cont_scale.sum(axis=0)
            mean = np.mean(summed_line[rms_start:rms_stop], dtype=np.float64)
            sigma = np.std(summed_line[rms_start:rms_stop], dtype=np.float64)
            rmsvalues[i] = sigma/mean
            print("Fringes scaled  {:10.3f}: RMS = {:8.4f}".format(current_scale[i], rmsvalues[i]))

        #
        # Determine scale factor that delivers the best RMS by an inverse-weighted average
        weight = 1.0/(rmsvalues*rmsvalues)
        rownum = range(nscales)
        weighted_rownum = weight*rownum
        minrms = rmsvalues.min()
        min_row = np.where(rmsvalues == minrms)[0][0]

        first_row = min_row - 2
        last_row = min_row + 3
        if min_row == 0 or min_row == nscales - 1:
            print(" ")
            print("WARNING: Best scale found on the edge of the specified scale range.")
            print("You are advised to try again after adjusting the scale range accordingly")
            thescale = current_scale[min_row]

        else:
            if min_row >= 2 and min_row <= nscales - 3:
                first_row = min_row - 2
                last_row = min_row + 3
            elif min_row == 1 or min_row == nscales - 2:
                first_row = min_row - 1
                last_row = min_row + 2
            w_shift_av = weighted_rownum[first_row:last_row].sum()
            weight_av = weight[first_row:last_row].sum()
            thescale = w_shift_av/weight_av
            thescale = thescale * scale_step + beg_scale

        print(" ")
        print(" Best scale : {:10.3f}".format(thescale))

        # Apply the best scale and create output array
        flat_scaled = fltdata.copy()
        flat_scaled[:, colstart+1:colstop-1] = (flat_scaled[:, colstart+1:colstop-1] - flat_mean) * thescale + \
                                               flat_mean

        # Write out file
        fits.writeto(outflat, data=flat_scaled)

        print("Output flat : {}".format(outflat))
        print("  (to be used as input to task 'defringe.py')")

    return

def get_flat_data(inflat, shifted_flat):
    if shifted_flat is not None:
        return shifted_flat
    else:
        f1 = fits.open(inflat)
        nextend = len(f1)
        if nextend == 0:
            return f1[0].data
        else:
            return f1[1].data


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
