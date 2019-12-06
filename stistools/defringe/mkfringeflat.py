from astropy.io import fits
from astropy.nddata.utils import block_reduce
import numpy as np
import math
import os
from scipy.ndimage import shift
import findloc
from response import response

__version__ = 0.1

def mkfringeflat(inspec, inflat, outflat, do_shift=True, beg_shift=-0.5, end_shift=0.5,
                 shift_step=0.1, do_scale=True, beg_scale=0.8, end_scale=1.2, scale_step=0.04,
                 extrloc=None, extrsize=None, opti_sreg=None, rms_region=None):

    print("mkfringeflat.py version {}".format(__version__))
    print(" - matching fringes in a flatfield to those in science data")

    sci_hdulist = fits.open(inspec)
    flt_hdulist = fits.open(inflat)

    scidata = sci_hdulist[1].data
    sci_prihdr = sci_hdulist[0].header
    sci_hdr = sci_hdulist[1].header

    if len(flt_hdulist) == 1:
        fltdata = flt_hdulist[0].data
        flthdr = flt_hdulist[0].header
    else:
        fltdata = flt_hdulist[1].data
        flthdr = flt_hdulist[1].header

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

    sax0 = round(-ltv1) + 1
    sax1 = sax0 + ncols - 1

    if opt_elem == "G750M":

        centera2 = sci_prihdr['centera2']
        sizaxis2 = sci_prihdr['sizaxis2']

        ltv2 = 1 - centera2 - sizaxis2/2
    
        flt_centera2 = flt_prihdr['centera2']
        flt_sizaxis2 = flt_flt_prihdr['sizaxis2']

        flt_ltv2 = 1 - flt_centera2 - flt_sizaxis2/2

    else:

        ltv2 = sci_hdr['ltv2']
        flt_ltv2 = flt_prihdr['sizaxis2']

    say0 = int(round(flt_ltv2 - ltv2))
    say1 = say0 + nrows - 1

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
        maxrow = int(round(findloc.find_loc(scidata)))
    else:
        maxrow = extrloc

    print(" Extraction center: row {}".format(maxrow))
    print("   Extraction size: {} pixels {}".format(numextrows, apername))
    fline = maxrow - int(round((numextrows - 0.49999)/2.0))
    lline = maxrow + int(round((numextrows - 0.49999)/2.0)) + 1

    if opti_sreg is None:
        if "G750M" in opt_elem:
            colstart = int(83/bincols)
            colstop = int(1106/bincols)
        else:
            colstart = int(5/bincols) - 1
            colstop = int(1020/bincols)

    if rms_region is None:
        rms_start = int(725/bincols) - 1
        rms_stop = int(900/bincols)

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
        nshifts = (lshift - fshift)/step + 1
        
        rmsvalues = np.zeros(nshifts)
        current_shift = np.zeros(nshifts)

        for i in range(nshifts):
            current_shift[i] = beg_shift + i*shift_step

            shifted_flat = shift(flt_blk, current_shift[i], order=1, mode='nearest')

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
            print("shift = {}, rms = {}".format(current_shift[i], rmsvalues[i]))

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
        print(" Best shift : {} pixels".format(theshift))

        # Apply the best shift and create output array
        flt_blk = block_reduce(fltdata, (binlines/fltbinlines, bincols/fltbincols),
                                func=np.mean)
        shifted_flat = shift(flt_blk, theshift, order=1, mode='nearest')

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
            print("Fringes scaled  {}: RMS = {}".format(current_scale[i], rmsvalues[i]))

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
            thescale = thescale*scale_step + beg_scale

        print(" ")
        print(" Best scale : {}".format(thescale))

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
