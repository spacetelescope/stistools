#! /usr/bin/env python

import sys
import os
import os.path
import math

import numpy as N
from scipy import signal as convolve

from astropy.io import fits
from . import gettable
from . import wavelen
from . import r_util

__version__ = "1.3 (2016 Feb 24)"


def wx2d(input, output, wavelengths=None, helcorr="",
         algorithm="wavelet",
         trace=None, order=7, subdiv=8, psf_width=0., rows=None,
         subsampled=None, convolved=None):
    """Resample the input, correcting for geometric distortion.

    Parameters
    ----------
    input : string
        name of input file containing an image set
    output : string
        name of the output file
    wavelengths : string, optional [Default: None]
        name of the output file for wavelengths
    helcorr : string
        specify "perform" or "omit" to override header keyword
    algorithm : {'wavelet', 'kd'}
        algorithm to use in resampling the input
    trace : string or array, or None
        trace array, or name of FITS table containing trace(s)
    order : int [Default: 7]
        polynomial order (an odd number, e.g. 5 or 7)
    subdiv : int [Default: 8]
        number of subpixels (a power of 2, e.g. 8 or 16)
    psf_width : float [Default: 0.]
        width of PSF for convolution (e.g. 1.3);
        0 means no convolution
    rows : tuple, optional [Default: None]
        a tuple giving the slice of rows to process; output values
        in all other rows will be set to zero.
        The default of None means all rows, same as (0, 1024)
    subsampled : string, optional [Default: None]
        name of the output file with the subsampled image
    convolved : string, optional [Default: None]
        name of the output file with the convolved image

    """

    if algorithm != "wavelet" and algorithm != "kd":
        raise ValueError("algorithm can only be 'wavelet' or 'kd'")
    if psf_width <= 0. and convolved is not None:
        print("Warning:  'convolved' will be ignored because psf_width=0")
        convolved = None
    if algorithm != "wavelet" \
       and (subsampled is not None or convolved is not None):
        raise ValueError(
            "cannot specify 'subsampled' or 'convolved' if algorithm='kd'")

    helcorr = helcorr.upper()

    ft = fits.open(input)
    nextend = ft[0].header.get("nextend", default=3)
    n_imsets = nextend // 3     # number of image sets

    if ft[1].header.get("ltm2_2", default=1.0) != 1.0:
        raise RuntimeError("can't handle data binned in the Y direction")

    phdu = ft[0]                # primary header/data unit

    # If trace was not specified, get the file name and selection keywords
    # from the header.
    tracefile = trace_name(trace, phdu.header)

    # Update the primary header, in preparation for writing to a new file.
    phdu.header.set("WX2DCORR", "COMPLETE", "this file is output from wx2d",
                    before="X2DCORR")
    is_an_array = isinstance(tracefile, N.ndarray)
    if is_an_array:
        phdu.header.add_history("trace array was specified explicitly")
    else:
        phdu.header.add_history("trace file " + tracefile)
    phdu.header.add_history("order = " + repr(order))
    phdu.header.add_history("subdiv = " + repr(subdiv))
    phdu.header.add_history("psf_width = " + repr(psf_width))

    # Write the primary header to the output file(s);
    # we'll append each extension in wx2d_imset.
    phdu.header.set("nextend", 0)           # no extensions yet
    phdu.header.set("filename", os.path.basename(output))
    phdu.writeto(output)
    if wavelengths is not None:
        phdu.header.set("filename", os.path.basename(wavelengths))
        phdu.writeto(wavelengths)
    if subsampled is not None:
        phdu.header.set("filename", os.path.basename(subsampled))
        phdu.writeto(subsampled)
    if convolved is not None:
        phdu.header.set("filename", os.path.basename(convolved))
        phdu.writeto(convolved)

    for imset0 in range(n_imsets):
        wx2d_imset(ft, imset0+1, output, wavelengths, helcorr,
                   algorithm,
                   tracefile, order, subdiv, psf_width, rows,
                   subsampled, convolved)

    ft.close()


def wx2d_imset(ft, imset, output, wavelengths, helcorr,
               algorithm,
               tracefile, order, subdiv, psf_width, rows,
               subsampled, convolved):
    """Resample one image set, and append to output file(s).

    Parameters
    ----------
    ft : HDUList
        Fits HDUList object for the input file.

    imset  :  int
        one-indexed image set number
    output  :  string
        name of the output file
    wavelengths :  string or None
        name of the output file for wavelengths
    helcorr :  {'perform', 'omit'}
        specify "perform" or "omit" to override header keyword
    algorithm : {"wavelet","kd"}
        algorithm to use to process input
    tracefile :  string or array
        trace array, or name of FITS table containing trace(s)
    order :  int
        polynomial order
    subdiv :  int
        number of subpixels
    psf_width :  float
        width of PSF for convolution
    rows :  tuple
        a tuple giving the slice of rows to process
    subsampled :  string, or None
        name of the output file with the subsampled image
    convolved :  string, or None
        name of the output file with the convolved image

    """

    hdu = ft[("SCI", imset)]
    errhdu = ft[("ERR", imset)]
    header = hdu.header
    nrows, ncols = hdu.data.shape
    original_nrows = nrows

    if rows is None:
        rows = (0, nrows)
    else:
        rows = (max(rows[0], 0), min(rows[1], nrows))
    nrows = rows[1] - rows[0]
    if original_nrows > nrows and imset == 1:
        ft[0].header.add_history("rows from %d to %d" % (rows[0]+1, rows[1]))

    # extract the subset
    if original_nrows > nrows:
        img = hdu.data[rows[0]:rows[1]]
        errimg = errhdu.data[rows[0]:rows[1]]
    else:
        img = hdu.data
        errimg = errhdu.data

    # offset of a row from nominal
    shifta2 = header.get("shifta2", default=0.)
    # offset to be applied to the trace array
    offset = rows[0] - header.get("ltv2", default=0.) + shifta2

    # Read the array of spectral traces (a2displ), and bin if the input
    # data are binned in the dispersion direction.
    (a2center, a2displ) = get_trace(tracefile, ft[0].header, header)

    if algorithm == "wavelet":
        (hdu.data, errhdu.data) = wavelet_resampling(hdu, img, errimg,
                                                     original_nrows, nrows, ncols, rows,
                                                     a2center, a2displ, offset, shifta2,
                                                     imset, order, subdiv, psf_width,
                                                     subsampled, convolved)
    else:
        (hdu.data, errhdu.data) = kd_resampling(img, errimg,
                                                original_nrows, nrows, ncols, rows,
                                                a2center, a2displ, offset, shifta2)
    del img

    # Write the SCI and ERR HDUs to the output file.
    ofd = fits.open(output, mode="update")
    ofd.append(hdu)
    ofd.append(errhdu)
    ofd.close()

    if wavelengths is not None:
        # Write an array of wavelengths.
        wl_hdu = fits.ImageHDU(header=hdu.header)
        if not helcorr:
            helcorr = ft[0].header.get("helcorr", default="OMIT").upper()
        if ft[0].header.get("sclamp", default="NONE") != "NONE":
            helcorr = "OMIT"
        wl_hdu.data = wavelen.compute_wavelengths((original_nrows, ncols),
                        ft[0].header, header, helcorr)
        ofd = fits.open(wavelengths, mode="update")
        ofd[0].header.set("nextend", imset)
        if helcorr == "PERFORM":
            ofd[0].header.set("helcorr", "COMPLETE")
        else:
            ofd[0].header.set("helcorr", "OMIT")
        ofd.append(wl_hdu)
        ofd.close()

    #  Create data quality array.

    hdu = ft[("DQ", imset)]
    im3 = hdu.data[rows[0]:rows[1]]
    nrows, ncols = im3.shape
    image3 = N.zeros((subdiv*nrows, ncols), dtype=im3.dtype)
    for j in range(subdiv):
        image3[j::subdiv, :] = im3
    del im3
    hdu.data[:, :] = 4           # "bad detector pixel or beyond aperture"
    hdu.data[rows[0]:rows[1]] = \
        apply_trace(image3, a2center, a2displ, subdiv, offset, shifta2, "DQ")
    del image3

    # Write the DQ HDU to the output file.
    ofd = fits.open(output, mode="update")
    ofd[0].header.set("nextend", imset*3)
    ofd.append(hdu)
    ofd.close()


def wavelet_resampling(hdu, img, errimg,
                       original_nrows, nrows, ncols, rows,
                       a2center, a2displ, offset, shifta2,
                       imset, order, subdiv, psf_width,
                       subsampled, convolved):

    """Resample img and errimg using wavelets.

    Parameters
    ----------
    hdu : fits header/data unit object
        header/data unit for a SCI extension
    img : ndarray
        SCI image array (could be a subset of full image)
    errimg : ndarray
        ERR image array (could be a subset of full image)
    original_nrows : int
        number of image lines (NAXIS2) in input image
    nrows : int
        number of image lines in subset
    ncols : int
        number of image columns (NAXIS1)
    rows : tuple
        tuple giving the slice of rows to process
    a2center : ndarray
        1-D array of Y locations
    a2displ : ndarray
        array of traces, one for each a2center; the length of each
        trace must be the same as the number of columns in the input image
    offset : float
        offset of the first row in 'image' from the beginning of
        the data block in the original file, needed for trace
    shifta2 : float
        offset of the row from nominal (from shifta2 keyword)
    imset : int
        number of the current image set (keyword EXTVER)
    order : int
        polynomial order
    subdiv : int
        number of subpixels per input pixel
    psf_width : float
        width of PSF for convolution (e.g. 1.3);
    subsampled : string or None
        name of the output file with the subsampled image
    convolved : string or None
        name of the output file with the convolved image

    Returns
    --------
    img_arr: tuple of ndarrays
        the image and error arrays (to replace the input img and errimg)

    """

    sub5 = N.zeros((nrows*subdiv, ncols), dtype=N.float32)
    err5 = N.zeros((nrows*subdiv, ncols), dtype=N.float32)
    for j in range(0, subdiv):
        sub5[j::subdiv] = img
        err5[j::subdiv] = errimg

    #  Wavelet interpolation loop
    # With each iteration, the number of rows is doubled and the values
    # in the tmp array (the subsampled SCI data) are halved; the SCI data
    # will later be extracted by summing subsampled values along a column.
    # In contrast, the ERR data will be extracted by taking the square root
    # of the average of the squares of the err array, so the err array
    # values should remain approximately constant rather than halving with
    # each iteration; this is the explanation for the extra factor of 2.

    # For CCD data the minimum value in the errimg array should be the
    # readout noise divided by the square root of the number of images
    # that have been combined, but for MAMA data it could be 0.
    # minerr is this minimum value, but set a lower limit of 1.
    minerr = max(1., errimg.min())

    step = subdiv
    while step > 1:

        tmp = N.zeros((nrows*2*subdiv//step, ncols), dtype=N.float32)
        err = tmp.copy()
        tmp[::2] = sub5[::step]
        tmp[1::2] = inv_avg_interp(order, tmp[::2])
        tmp = inv_haar(tmp)
        # for computing the error, use minerr as the minimum flux
        err_tmp0 = N.maximum(tmp[0::2], minerr)
        err_tmp1 = N.maximum(tmp[1::2], minerr)
        err[0::2] = N.sqrt(2.*err_tmp0 / (err_tmp0 + err_tmp1)) * err5[::step]
        err[1::2] = N.sqrt(2.*err_tmp1 / (err_tmp0 + err_tmp1)) * err5[::step]
        step //= 2

        sub5[::step] = tmp
        err5[::step] = err
        for j in range(step):
            sub5[j::step] = tmp
            err5[j::step] = err

    del (err_tmp0, err_tmp1)

    if subsampled is not None:
        hdu.data = sub5.copy()
        ofd = fits.open(subsampled, mode="update")
        ofd[0].header.set("nextend", imset)
        ofd.append(hdu)
        ofd.close()

    #  Optional PSF convolution

    if psf_width > 0.:
        cnv = N.zeros(sub5.shape, dtype=N.float32)
        krn = N.array([stis_psf(float(j), psf_width*subdiv)
                       for j in range(-32, 32)])
        krn = krn / N.sum(krn)
        for j in range(sub5.shape[1]):
            cnv[:, j] = convolve.convolve(sub5[:, j], krn, mode='same')  # same size
        if convolved is not None:
            hdu.data = cnv.copy()
            ofd = fits.open(convolved, mode="update")
            ofd[0].header.set("nextend", imset)
            ofd.append(hdu)
            ofd.close()
    else:
        cnv = sub5

    if original_nrows > nrows:
        result = N.zeros((original_nrows, ncols), dtype=N.float32)
        err_result = N.zeros((original_nrows, ncols), dtype=N.float32)
        result[rows[0]:rows[1]] = apply_trace(cnv, a2center, a2displ,
                                              subdiv, offset, shifta2, "SCI")
        err_result[rows[0]:rows[1]] = apply_trace(err5, a2center, a2displ,
                                                  subdiv, offset, shifta2, "ERR")
    else:
        result = apply_trace(cnv, a2center, a2displ,
                             subdiv, offset, shifta2, "SCI")
        err_result = apply_trace(err5, a2center, a2displ,
                                 subdiv, offset, shifta2, "ERR")

    return result, err_result


def kd_resampling(img, errimg,
                  original_nrows, nrows, ncols, rows,
                  a2center, a2displ, offset, shifta2):
    """Apply Kris Davidson's resampling method.

    Parameters
    ----------
    img :  ndarray
        SCI image array (could be a subset of full image)
    errimg :  ndarray
        ERR image array (could be a subset of full image)
    original_nrows :  int
        number of image lines (NAXIS2) in input image
    nrows :  int
        number of image lines in subset
    ncols :  int
        number of image columns (NAXIS1)
    rows :  tuple
        tuple giving the slice of rows to process
    a2center : ndarray
        1-D array of Y locations
    a2displ :  ndarray
        array of traces, one for each a2center; the length of each
        trace must be the same as the number of columns in the input image
    offset : float
        offset of the first row in 'image' from the beginning of
        the data block in the original file, needed for trace
    shifta2 : float
        offset of the row from nominal (from shifta2 keyword)

    Returns
    -------
    img_arr : tuple
        the image and error arrays (to replace the input img and errimg)

    """

    # image2 is for the error image
    subdiv = 8
    image2 = N.zeros((subdiv*nrows, ncols), dtype=N.float32)
    for j in range(subdiv):
        image2[j::subdiv, :] = errimg

    if original_nrows > nrows:
        result = N.zeros((original_nrows, ncols), dtype=N.float32)
        err_result = N.zeros((original_nrows, ncols), dtype=N.float32)
        result[rows[0]:rows[1]] = kd_apply_trace(img, a2center, a2displ,
                                                 offset, shifta2)
        err_result[rows[0]:rows[1]] = apply_trace(image2, a2center, a2displ,
                                                  subdiv, offset, shifta2, "ERR")
    else:
        result = kd_apply_trace(img, a2center, a2displ, offset, shifta2)
        err_result = apply_trace(image2, a2center, a2displ,
                                 subdiv, offset, shifta2, "ERR")

    return result, err_result


def kd_apply_trace(image, a2center, a2displ, offset=0., shifta2=0.):
    """Kris Davidson's resampling algorithm, following the trace.

    Parameters
    ----------
    image :  ndarray
        input 2-D image array
    a2center : ndarray
        array of Y locations
    a2displ : ndarray
        array of traces, one for each a2center; the length of
        each trace must be the same as the number of columns in 'image'
    offset : float
        offset of the first row in 'image' from the beginning
        of the data block in the original file, needed for trace
    shifta2 : float
        offset of the row from nominal (from shifta2 keyword)

    Returns
    -------
    x2d : ndarray
        2-D array containing the resampled image

    """

    shape = image.shape
    x2d_shape = N.array(shape)
    x2d = N.zeros(x2d_shape, dtype=N.float32)
    total = shape[0] * shape[1]

    flat_im = N.ravel(image)

    for i in range(x2d_shape[0]):
        # y is the location in the cross-dispersion direction.
        trace = interpolate_trace(a2center, a2displ, float(i) + offset,
                                  x2d_shape[1])
        y = float(i) + shifta2 + trace

        nint_y = N.around(y)
        s = y - nint_y
        n = nint_y.astype(N.int32)
        col_range = N.arange(x2d_shape[1], dtype=N.int32)

        # These are indices into the flattened image.
        nm2 = (n-2) * x2d_shape[1] + col_range
        nm1 = (n-1) * x2d_shape[1] + col_range
        n0 = n * x2d_shape[1] + col_range
        np1 = (n+1) * x2d_shape[1] + col_range
        np2 = (n+2) * x2d_shape[1] + col_range

        nm2 = N.where(nm2 < 0, 0, nm2)
        nm1 = N.where(nm1 < 0, 0, nm1)
        n0 = N.where(n0 < 0, 0, n0)
        np1 = N.where(np1 < 0, 0, np1)
        np2 = N.where(np2 < 0, 0, np2)

        nm2 = N.where(nm2 >= total, total-1, nm2)
        nm1 = N.where(nm1 >= total, total-1, nm1)
        n0 = N.where(n0 >= total, total-1, n0)
        np1 = N.where(np1 >= total, total-1, np1)
        np2 = N.where(np2 >= total, total-1, np2)

        a = - 0.050 * flat_im[nm2] + 0.165 * flat_im[nm1] + 0.77 * flat_im[n0] \
            + 0.165 * flat_im[np1] - 0.050 * flat_im[np2]
        b = 0.08 * flat_im[nm2] - 0.66 * flat_im[nm1] \
            + 0.66 * flat_im[np1] - 0.08 * flat_im[np2]
        c = 0.04 * flat_im[nm2] + 0.34 * flat_im[nm1] - 0.76 * flat_im[n0] \
            + 0.34 * flat_im[np1] + 0.04 * flat_im[np2]

        x2d[i] = a + b * s + c * s**2

    return x2d


def stis_psf(x, a):

    """Evaluate the cross-dispersion PSF at x.

    Parameters
    ----------
    x : float
        offset in pixels from the center of the profile
    a : float
        a measure of the width of the PSF

    Returns
    -------
    val : float
        the PSF evaluated at x

    """

    return (1. + (x/float(a))**2)**-2


def apply_trace(image, a2center, a2displ, subdiv,
                offset=0., shifta2=0., extname="SCI"):
    """Add together 'subdiv' rows of 'image', following the trace.

    Parameters
    ----------
    image :  ndarray
        input 2-D image array, oversampled by 'subdiv' in axis 0
    a2center : ndarray
        1-D array of Y locations
    a2displ : ndarray
        array of traces, one for each a2center; the length of each
        trace must be the same as the number of columns in the input image
    subdiv :  int
        number of rows to add together
    offset :  float
        offset of the first row in 'image' from the beginning of
        the data block in the original file, needed for trace
    shifta2 :  float
        offset of the row from nominal (from shifta2 keyword)
    extname :  string
        which type of extension (SCI, ERR, DQ)?

    Returns
    ---------
    x2d : ndarray
        resampled 2-D image array

    Notes
    -------
    The function value is a 2-D array containing the resampled image.
    This is binned by subdiv in Y (axis 0), after shifting by trace
    (multiplied by subdiv).

    For extname = "ERR" the result differs in these ways:

      (1) fractions of pixels at the endpoints of the extraction region
          are not included
      (2) the values are combined as the average of the sum of the squares

    For extname = "DQ" the result differs in these ways:

      (1) the output is type int16
      (2) the output values are nominally the same as the input, while
          for SCI the output are subdiv times larger than the input
      (3) fractions of pixels at the endpoints of the extraction region
          are not included
      (4) the values are combined via bitwise OR rather than an average or sum

    """

    shape = image.shape
    x2d_shape = N.array(shape)
    x2d_shape[0] //= subdiv
    if extname == "DQ":
        x2d = N.zeros(x2d_shape, dtype=N.int16)
    else:
        x2d = N.zeros(x2d_shape, dtype=N.float32)

    for i in range(x2d_shape[0]):
        # y is the location in the output, binned image (x2d), while
        # locn is the location in the input, oversampled image.
        trace = interpolate_trace(a2center, a2displ, float(i) + offset,
                                  x2d_shape[1])
        y = float(i) + shifta2 + trace
        locn = y * subdiv + (subdiv - 1.) / 2.
        if extname == "SCI":
            x2d[i] = extract(image, locn, subdiv)
        elif extname == "ERR":
            x2d[i] = extract_err(image, locn, subdiv)
        else:
            x2d[i] = extract_i16(image, locn, subdiv)

    return x2d


def extract(image, locn, subdiv):
    """Add together 'subdiv' rows of 'image', centered on 'locn'.

    Parameters
    ----------
    image : ndarray
        input array, oversampled by 'subdiv' in axis 0
    locn : ndarray
        a 1-D array giving the location at which to extract; an
        integer value corresponds to the center of the pixel.  The length
        must be the same as the number of columns in the input image.
    subdiv : int
        number of rows to add together

    Returns
    --------
    spec : ndarray
        a 1-D array containing the extracted row

    """

    # Shift the zero point so the edges of a pixel have integer coordinates.
    locn0 = locn + 0.5

    shape = image.shape

    hw = subdiv // 2
    fhw = float(hw)
    # floating point locations of upper and lower edges
    fhigh = locn0 + fhw
    flow = locn0 - fhw
    # integer endpoints of range of whole pixels
    high = N.floor(fhigh).astype(N.int32)
    low = N.ceil(flow).astype(N.int32)
    # fractions of pixels at upper and lower edges
    dhigh = fhigh - high
    dlow = low - flow

    spec = N.zeros(shape[1], dtype=N.float32)
    for j in range(shape[1]):
        s_low = low.item(j)
        s_high = high.item(j)
        if s_low < 1 or s_high > shape[0]-1:
            continue
        spec[j] = N.sum(image[s_low: s_high, j])
        spec[j] += image.item(s_high, j) * dhigh.item(j)
        spec[j] += image.item(s_low-1, j) * dlow.item(j)

    return spec


def extract_err(image, locn, subdiv):
    """Average 'subdiv' rows of 'image', centered on 'locn'.

    Parameters
    ----------
    image : ndarray
        input array, oversampled by 'subdiv' in axis 0
    locn : ndarray
        a 1-D array giving the location at which to extract; an
        integer value corresponds to the center of the pixel
    subdiv : int
        number of rows to add together

    Returns
    -------
    spec : ndarray
        a 1-D array containing the extracted row

    Notes
    -------
    This takes the square root of the average of the squares, intended to be
    used for interpolating the ERR array.  Fractions of pixels at the upper
    and lower edges are excluded.
    """

    # Shift the zero point so the edges of a pixel have integer coordinates.
    locn0 = locn + 0.5

    shape = image.shape

    hw = subdiv // 2
    fhw = float(hw)
    # floating point locations of upper and lower edges
    fhigh = locn0 + fhw
    flow = locn0 - fhw
    # integer endpoints of range of whole pixels
    high = N.floor(fhigh).astype(N.int32)
    low = N.ceil(flow).astype(N.int32)

    spec = N.zeros(shape[1], dtype=N.float32)
    for j in range(shape[1]):
        s_low = low.item(j)
        s_high = high.item(j)
        if s_low < 1 or s_high > shape[0]-1:
            continue
        sum = 0.
        for i in range(s_low, s_high+1):
            s_image = image.item(i, j)
            sum += s_image**2
        sum /= (s_high - s_low + 1.)
        spec[j] = math.sqrt(sum)

    return spec


def extract_i16(image, locn, subdiv):
    """Bitwise OR 'subdiv' rows of 'image', centered on 'locn'.

    Parameters
    ----------
    image : ndarray
        input array, oversampled by 'subdiv' in axis 0
    locn : ndarray
        a 1-D array giving the location at which to extract; an
        integer value corresponds to the center of the pixel
    subdiv : int
        number of rows to add together

    Returns
    --------
    spec : ndarray
        a 1-D array containing the extracted row

    """

    # Shift the zero point so the edges of a pixel have integer coordinates.
    locn0 = locn + 0.5

    shape = image.shape

    hw = subdiv // 2
    fhw = float(hw)
    # floating point locations of upper and lower edges
    fhigh = locn0 + fhw
    flow = locn0 - fhw
    # integer endpoints of range of whole pixels
    high = N.floor(fhigh).astype(N.int32)
    low = N.ceil(flow).astype(N.int32)

    spec = N.zeros(shape[1], dtype=N.int16)
    for j in range(shape[1]):
        s_low = low.item(j)
        s_high = high.item(j)
        if s_low < 1 or s_high > shape[0]-1:
            continue
        sum = 0
        for i in range(s_low, s_high+1):
            sum |= image.item(i, j)
        spec[j] = sum

    return spec


def interpolate_trace(a2center, a2displ, y, length):
    """Interpolate within the array of traces, and return a trace.

    Parameters
    ----------
    a2center : ndarray
        array of Y locations
    a2displ : ndarray
        array of traces, one trace for each element of a2center
    y : float
        Y location on the detector
    length : int
        length of a trace; needed only if traces is empty

    """

    if len(a2displ) < 1:
        trace = N.zeros(length, dtype=N.float32)
    else:
        # interpolate to get the trace at y
        trace = r_util.interpolate(a2center, a2displ, y)

    return trace


def trace_name(trace, phdr):
    """Return the 1dt table name or array.

    Parameters
    -----------
    trace : string or array or None
        if trace is None the header keyword SPTRCTAB will be
        gotten from phdr; else if this is a string it should be the name
        of a trace file (possibly using an environment variable); otherwise,
        it should be a trace, in which case it will be returned unchanged
    phdr : fits Header object
        primary header, used only if trace is None

    Returns
    --------
    tracefile : string or array
        name of a trace file (with environment variable expanded),
        or an actual trace array

    """

    if trace is None:
        try:
            tracefile = phdr["sptrctab"]
        except KeyError:
            raise ValueError("Keyword SPTRCTAB not found; specify trace explicitly.")
        tracefile = r_util.expandFileName(tracefile)
    elif isinstance(trace, str):
        tracefile = r_util.expandFileName(trace)
    else:
        # trace could already be an array object
        tracefile = trace

    return tracefile


def get_trace(tracefile, phdr, hdr):
    """Read 1-D traces from the 1dt table (sptrctab).

    Parameters
    ----------
    tracefile : string or array
        either a trace array or the name of a FITS 1dt table
    phdr : fits Header object
        primary header of input file
    hdr : fits Header object
        extension header of input image (for binning info and
        time of exposure)

    Returns
    --------
    trace_arrays : tuple of 2 arrays
        a pair of arrays, one is the Y location at the middle column,
        and the other is an array of trace arrays

    Notes
    -----
    If 'tracefile' is already a trace array, it will just be returned,
    together with an arbitrary Y location of 0 (because that will
    always be within the image).

    opt_elem and cenwave are criteria for selecting the relevant rows
    from the 1dt table.  There will normally be several rows that match,
    and they should have different values of the Y location; the output
    list will be sorted on Y location.
    """

    # If the input is already an array object, return a2center
    # (arbitrarily set to 0) and the trace array.
    is_an_array = isinstance(tracefile, N.ndarray)
    if is_an_array:
        a2center = [0.]
        a2displ = [tracefile]
        return a2center, a2displ

    # These are the row-selection values.
    opt_elem = phdr.get("opt_elem", default="")
    cenwave = phdr.get("cenwave", default=0)

    # These describe the binning in the dispersion direction.
    ltm = hdr.get("ltm1_1", default=1.)
    ltv = hdr.get("ltv1", default=0.)                  # one indexing
    binaxis1 = int(round(1. / ltm))
    if binaxis1 not in [1, 2, 4, 8]:
        raise ValueError("LTM1_1 should be either 1., 0.5, 0.25, or 0.125")

    filter = {"opt_elem": opt_elem, "cenwave": cenwave}
    trace_info = gettable.getTable(tracefile, filter,
                                   sortcol="a2center", at_least_one=True)
    expstart = hdr.get("expstart", default=-1.)
    gettable.rotateTrace(trace_info, expstart)
    a2center = trace_info.field("a2center") - 1.       # zero indexing
    a2displ = trace_info.field("a2displ")

    if binaxis1 > 1 and len(a2displ) > 0:
        a2displ = bin_traces(a2displ, binaxis1, ltv)

    return a2center, a2displ


def bin_traces(a2displ, binaxis1, ltv):
    """bin the traces by the factor binaxis1

    Parameters
    ----------
    a2displ : ndarray
        an array of one or more arrays of Y displacements (traces)
    binaxis1 : int
        binning factor in the dispersion axis
    ltv : float
        offset in the dispersion axis (one indexing)

    Returns
    -------
    a2displ : ndarray
        an array of traces (a2displ), but with the trace arrays binned
        and shorter by the factor binaxis1

    """

    if len(a2displ) < 1 or binaxis1 == 1:
        return a2displ

    oldlen = len(a2displ[0])
    if binaxis1 == 2:
        newlen = 511
    elif binaxis1 == 4:
        newlen = 255
    elif binaxis1 == 8:
        newlen = 127
    else:
        newlen = oldlen
    newtraces = []

    # Note:  one-indexed pixels are used for determining the offset.
    # The left edge of the first image pixel is pixel number 0.5; convert
    # that point in the binned image to reference pixel coordinates.
    # Then compare with the left edge of the first reference pixel to get
    # the offset of our binned image (in reference pixel units).
    left_edge = (0.5 - ltv) * binaxis1
    offset = left_edge - 0.5                    # should be an integer
    offset = int(round(offset))

    for trace in a2displ:
        newtrace = N.zeros(newlen, dtype=N.float32)
        for i in range(binaxis1):
            # this slice can be longer than newlen
            temp = trace[i+offset:oldlen:binaxis1]
            newtrace += temp[0:newlen]
        newtrace /= float(binaxis1)
        newtraces.append(newtrace)

    return N.array(newtraces)


def inv_haar(image):
    image[0::2] = (image[0::2] + image[1::2]) / 2.
    image[1::2] = (image[0::2] - image[1::2])
    return image


def inv_avg_interp(order, image):

    side, j0, j1 = (order-1)//2, (order-1)//2, (order+1)//2
    rows, cols = image.shape

    order_2 = float(order) / 2.
    n = order + 1

    x = N.arange(n, dtype=N.float64)
    y = N.zeros((n, cols), dtype=N.float64)
    d = N.zeros((rows, cols), dtype=N.float32)  # single precision result

    for j in range(side, rows-side):
        y[1:] = N.cumsum(image[j-side:j+side+1], axis=0)
        d[j] = -((y[j1] + y[j0]) - 2. * polynomial(x, y, order_2, n))
    return d


def polynomial(x, y, z, n):
    """used for interpolation

    Parameters
    ----------
    x : ndarray
        the integer values from 0 through n-1 inclusive (but float64)
    y : ndarray
        a 2-D array, axis 0 of length n
    z : float
        n / 2.
    n : int
        1 + order of polynomial fit
    """

    t = n*[0.]
    for k in range(n):
        t[k] = y[k]
        for j in range(k-1, -1, -1):
            t[j] = t[j+1] + (t[j+1] - t[j]) / \
                            ((z - x.item(j)) / (z - x.item(k)) - 1.)
    return t[0]

if __name__ == "__main__":

    # Note that the command-line options do not include all of the
    # wx2d function arguments.
    if len(sys.argv) < 3 or len(sys.argv) > 7:
        print("Syntax:  wx2d.py input output "
              "[wavelengths trace minrow maxrow]")
        sys.exit()

    input = sys.argv[1]
    output = sys.argv[2]
    wavelengths = None
    trace = None

    if len(sys.argv) > 3:
        if sys.argv[3] != "None":
            wavelengths = sys.argv[3]

    if len(sys.argv) > 4:
        if sys.argv[4] != "None":
            trace = sys.argv[4]

    if len(sys.argv) == 7:
        row0 = int(sys.argv[5]) - 1
        row1 = int(sys.argv[6])
        wx2d(sys.argv[1], sys.argv[2], wavelengths=wavelengths, trace=trace,
             rows=(row0, row1))
    else:
        wx2d(sys.argv[1], sys.argv[2], wavelengths=wavelengths, trace=trace)
