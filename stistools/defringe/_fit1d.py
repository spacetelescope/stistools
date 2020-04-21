#
# This module reproduces much of the behaviour of IRAF 1-d fitting routines

import math
import numpy as np
from scipy.interpolate import LSQUnivariateSpline

def fit1d(x, y, weights=None,
          naverage=1,
          function="spline3",
          order=3,
          low_reject=3.0,
          high_reject=3.0,
          niterate=0,
          grow=0.0):
    """Main 1-d fitting function.  The parameters can be used exactly like they are in IRAF tasks
    like fit1d, response, etc.  The default behaviour is to fit a cubic spline with order 3
    (3 equal regions, 2 knots) with no sigma clipping.

    Parameters
    ----------

    x: 1-d ndarray
        The array of independent values

    y: 1-d ndarray
        The array of dependent values.  Must have same length as x, or function returns None

    weights: 1-d ndarray or None
        The array of weights.  If set to None, the points are given equal weight of 1.0.  If not
        set to None, must have the same length as x and y, or the function returns None

    naverage: int
        The number of adjacent elements of x, y, and weight that are averaged before
        fitting is done.  If naverage=1, no averaging is done.

    function: str
        Fitting function.  Currently only "spline1" and "spline3" are supported

    order: int
        Order of fitting fuction.  The number of knots in the spline fit is (order-1)

    low_reject: float
        Points with y values less than fit - low_reject*(rms deviation of data from fit) are rejected
        when the fit is iterated.  Ignored if niterate = 0.

    high_reject: float
        Points with y values greater than fit + low_reject*(rms deviation of data from fit) are rejected
        when the fit is iterated.  Ignored if niterate = 0.

    niterate: int
        Number of times sigma-clipping iterations are performed.  niterate=0 corresponds to no
        sigma clipping

    grow: float
        Used to calculate how many elements adjacent to reject pixels are also rejected.
        Not currently implemented, but included to duplicate the IRAF fitting parameter.

    Returns
    -------

    fit: scipy.interpolate.fitpack2.LSQUnivariateSpline object
        The result of the fit.  It can be used to calculate the fitted
        values for the input x values:

        fitted_values = fit(x)

    """

    nx = len(x)
    ny = len(y)
    if nx != ny:
        print("X and Y vectors have different lengths, aborting fit1d")
        return None
    if weights is not None:
        nw = len(weights)
        if nw != nx:
            print("Weight vector has different length to data vectors, aborting fit1d")
            print("Use weight=None to assign equal weights")
            return None
    else:
        weights = np.ones(nx)
    xfit, yfit, wfit = wtrebin(x, y, weights, naverage)
    fitted = fit_with_rejection(xfit, yfit, wfit, function, order, low_reject, high_reject,
                                niterate, grow)
    return fitted

def get_knots(x, number_of_knots):
    """Calculate equally-spaced knots for an array, given the array and the number of knots
    to calculate.  The number of knots is 1 fewer than the order of the fit.

    Parameters
    ----------

    x: 1-d ndarray
        The input array for which knots are to be calculated

    number_of_knots: int
        The number of knots to calculate

    Returns
    -------

    knots: 1-d ndarray
        The array of equally-spaced knots

    """

    interval = x[-1] - x[0]
    #
    # If there are n knots, we are dividing the interval into (n+1) regions
    #
    subinterval = interval/float(number_of_knots + 1)
    knots = x[0] + subinterval*np.arange(number_of_knots + 1)
    #
    # First knot is at x[0], we don't need that one
    return knots[1:]

def fit_once(x, y, weights, function, order):
    """Do a single fit of function to (x, y) data.

    Parameters
    ----------

    x: 1-d ndarray
        Input array of x-values

    y: 1-d array
        Input array of y-values

    weights: 1-d array
        Input array of weights

    function: str
        Function to be fitted.  Currently only "spline1" and "spline3" are supported

    order: int
        Order of function to be fitted.  Since only splines are currently supported,
        the order is the same as the number of equal reqions the data domain is split into,
        and 1 more than the number of equally-spaced knots.

    Returns
    -------

    fit: scipy.interpolate.fitpack2.LSQUnivariateSpline object
        This can be used to calculate the fitted values for the input x values:

        fitted_values = fit(x)

    """

    nx = len(x)
    if function[:6] == "spline":
        #
        # Number of knots is 1 fewer than the order
        knots = get_knots(x, order-1)
        if function == "spline3":
            k = 3
        elif function == "spline1":
            k = 1
        else:
            print("Spline fitting only valid for spline1 and spline3 functions")
            return None
        fitted = LSQUnivariateSpline(x, y, knots, w=weights, k=k)
    else:
        print("Not implemented yet")
        return None
    return fitted

def fit_with_rejection(x, y, weights, function, order, low_reject, high_reject, niterate, grow):
    """Fit with sigma-clipping

    Parameters
    ----------

    x: 1-d ndarray
        The array of independent values

    y: 1-d ndarray
        The array of dependent values.  Must have same length as x, or function returns None

    weights: 1-d ndarray or None
        The array of weights.  If set to None, the points are given equal weight of 1.0.  If not
        set to None, must have the same length as x and y, or the function returns None

    function: str
        Fitting function.  Currently only "spline1" and "spline3" are supported

    order: int
        Order of fitting fuction.  The number of knots in the spline fit is (order-1)

    low_reject: float
        Points with y values less than fit - low_reject*(rms deviation of data from fit) are rejected
        when the fit is iterated.  Ignored if niterate = 0.

    high_reject: float
        Points with y values greater than fit + low_reject*(rms deviation of data from fit) are rejected
        when the fit is iterated.  Ignored if niterate = 0.

    niterate: int
        Number of times sigma-clipping iterations are performed.  niterate=0 corresponds to no
        sigma clipping

    grow: float
        Used to calculate how many elements adjacent to reject pixels are also rejected.
        Not currently implemented, but included to duplicate the IRAF fitting parameter.

    Returns
    -------

    fit: scipy.interpolate.fitpack2.LSQUnivariateSpline object
        The result of the fit.  It can be used to calculate the fitted
        values for the input x values:

        fitted_values = fit(x)

    """

    fitted = fit_once(x, y, weights, function, order)
    if niterate <= 0:
        return fitted
    npts = len(x)
    rej_weights = np.ones(npts)
    for iter in range(niterate):
        rms_deviation = calc_rms_deviation(x, y, weights*rej_weights, fitted)
        deviation = (y - fitted(x))*rej_weights
        above_hi_rej = np.where(deviation > high_reject*rms_deviation)
        rej_weights[above_hi_rej] = 0.0
        nhigh = len(above_hi_rej[0])
        below_low_rej = np.where(deviation < -low_reject*rms_deviation)
        rej_weights[below_low_rej] = 0.0
        nlow = len(below_low_rej[0])
        nbad = nhigh + nlow
        if nbad == 0: break
        fitted = fit_once(x, y, weights*rej_weights, function, order)
    return fitted

def calc_rms_deviation(x, y, weights, fitted):
    """Calculate the weighted RMS deviation between y and fitted(x)

    Parameters
    ----------

    x: 1-d ndarray
        Input array of x-values

    y: 1-d array
        Input array of y-values

    weights: 1-d array
        Input array of weights

    fitted: scipy.interpolate.fitpack2.LSQUnivariateSpline object
        The result from running the fitting routine

    Returns
    -------

    rms_deviation: float
        The rms deviation

    """

    sumsq = 0.0
    npts = len(x)
    deviation = y - fitted(x)
    devsquared = weights*deviation*deviation
    sumsq = devsquared.sum()
    sumweights = weights.sum()
    rms_deviation = (np.sqrt(sumsq/sumweights))
    return rms_deviation

def wtrebin(x, y, weights=None, nbin=1):
    """Like the IRAF rg_wtbin.  Bin the x, y, and weight data by averaging nbin values.

    Parameters
    ----------

    x: 1-d ndarray
        The array of independent values

    y: 1-d ndarray
        The array of dependent values.  Must have same length as x, or function
        returns (None, None, None)

    weights: 1-d ndarray or None
        The array of weights.  If set to None, the points are given equal weight of 1.0.  If not
        set to None, must have the same length as x and y, or the function
        returns (None, None, None)

    nbin: int
        The number of adjacent elements of x, y, and weight that are averaged before
        fitting is done.  If nbin=1, no averaging is done and the input arrays are returned

    Returns
    -------

    (xout, yout, wout): 3-tuple of ndarray objects
        The binned/averaged x, y, and weight arrays

    """

    if nbin == 1:
        return x, y, weights
    nx = len(x)
    ny = len(y)
    if nx != ny:
        print("Lengths of x and y vectors different, aborting")
        return None, None, None
    if type(weights) == type(None):
        weights = np.ones(nx)
    nw = len(weights)
    if nx != nw:
        print("Lengths of data and weight vectors different, aborting")
        return None, None, None
    nout = int(math.ceil(nx/float(nbin)))
    xout_temp = np.zeros((nbin, nout))
    yout_temp = np.zeros((nbin, nout))
    wtout_temp = np.zeros((nbin, nout))
    for i in range(nbin):
        ntemp = len(x[i::nbin])
        xout_temp[i, :ntemp] = x[i::nbin]*weights[i::nbin]
        yout_temp[i, :ntemp] = y[i::nbin]*weights[i::nbin]
        wtout_temp[i, :ntemp] = weights[i::nbin]
    xout = xout_temp.sum(axis=0)/wtout_temp.sum(axis=0)
    yout = yout_temp.sum(axis=0)/wtout_temp.sum(axis=0)
    wtout = wtout_temp.sum(axis=0)/float(nbin)
    return xout, yout, wtout
