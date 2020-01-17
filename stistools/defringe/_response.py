import math
import numpy as np
from astropy.io import fits

from ._fit1d import fit1d

def response(calibration_array, normalization_array,
             threshold=None,
             function="spline3",
             sample='*',
             naverage=2,
             order=10,
             low_reject=3.0,
             high_reject=3.0,
             niterate=2,
             grow=0.0):
    """
    RESPONSE -- Determine the response function for 2D spectra.
    
    A calibration image is divided by a normalization spectrum to form
    a response image.  The normalization spectrum is derived by averaging
    the normalization image across dispersion.  The normalization spectrum
    is then smoothed by curve fitting.  The smoothed normalization
    spectrum is divided into the calibration image to form the response
    function image.  The curve fitting may be performed interactively
    using the icfit package.  A response function is determined for each
    input image.  Image sections in the calibration image may be used to determine
    the response for only part of an image such as with multiple slits.
    
    The images are given by image templates.  The number of images must
    in each list must match.  Image sections are allowed in the calibration
    image.
    

    lpar response:
  calibration = "star_cont_sh1.fits[5:1020,505:515]" Longslit calibration images
 normalizatio = "star_cont_sh1.fits[5:1020,505:515]" Normalization spectrum images
     response = "norm_star_cont_sh1.fits" Response function images
 (interactive = yes)            Fit normalization spectrum interactively?
   (threshold = INDEF)          Response threshold
      (sample = "*")            Sample of points to use in fit
    (naverage = 2)              Number of points in sample averaging
    (function = "spline3")      Fitting function
       (order = 15)             Order of fitting function
  (low_reject = 3.0)            Low rejection in sigma of fit
 (high_reject = 3.0)            High rejection in sigma of fit
    (niterate = 2)              Number of rejection iterations
        (grow = 0.0)            Rejection growing radius
    (graphics = "stdgraph")     Graphics output device
      (cursor = "")             Graphics cursor input
        (mode = "al")           

    """

    if sample != "*":
        print ('Only sample="*" currently supported for this version of response')
    response = make_response(calibration_array, normalization_array, threshold, naverage,
                             function, order, low_reject, high_reject, niterate, grow)
    return response

def make_response(calibration_array, normalization_array, threshold, naverage,
                  function, order, low_reject, high_reject, niterate, grow):
    """Like re_make in response.x

    """

    dispaxis = get_dispaxis(calibration_array)

    # Get the normalization spectrum
    nrows, ncols = normalization_array.shape
    average_spectrum = normalization_array.sum(axis=0) / float(nrows)
    x = np.arange(ncols) + 1
    fitted = fit1d(x, average_spectrum, weights=None, function=function, order=order,
                   naverage=naverage, low_reject=low_reject, high_reject=high_reject,
                   niterate=niterate, grow=grow)

    fitted_spectrum = fitted(x)
    
    response = normalise_response(calibration_array, dispaxis, threshold, fitted_spectrum)

    return response

def get_dispaxis(calibration_array):
    """For completeness, include a function to get this.  For our purposes, dispaxis is
    always 1

    """

    return 1

def normalise_response(calibration_array, dispaxis, threshold, fitted_spectrum):
    """Like re_normalize in response.x

    """

    nrows, ncols = calibration_array.shape

    response = calibration_array / fitted_spectrum

    if threshold is not None:
        low_spectrum = np.where(fitted_spectrum < threshold)
        for i in low_spectrum[0]:
            response[:,i] = 1.0
        low_cal_array = np.where(calibration_array < threshold)
        response[low_cal_array] = 1.0

    return response
