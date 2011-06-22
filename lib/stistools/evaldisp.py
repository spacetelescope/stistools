from __future__ import division         # confidence high

def newton (x, coeff, cenwave, niter=4):
    """Return the wavelength corresponding to pixel x.

    The dispersion solution is evaluated iteratively, and the slope
    (dispersion) for Newton's method is determined numerically, using a
    difference in wavelength of one Angstrom.  Note that the evalDisp
    in this file assumes that the grating is first order.

    @param x:  a single pixel number or an array of pixel numbers
    @type x:  float or array of Float64

    @param coeff:  a list of eight elements containing the
        dispersion coefficients as read from a STIS _dsp.fits table
    @type coeff:  either a list of float or a numarray of Float64

    @param cenwave:  central wavelength, in Angstroms
    @type cenwave:  int or float

    @param niter:  number of iterations
    @type niter:  int

    @return:  a single wavelength or an array (numarray) of wavelengths,
        in Angstroms
    @rtype:  either a float or a numarray of Float64
    """

    wl = cenwave
    x0 = evalDisp (coeff, wl)

    delta_wl = 1.                       # one Angstrom
    for i in range (niter):
        x1 = evalDisp (coeff, wl+delta_wl)
        dispersion = delta_wl / (x1 - x0)
        wl += dispersion * (x - x0)
        x0 = evalDisp (coeff, wl)

    return wl

def evalDisp (coeff, wl):
    """Return the pixel corresponding to wavelength wl.

    The expression in the calstis code is:

    x = coeff[0] +
        coeff[1] * m * wl +
        coeff[2] * m**2 * wl**2 +
        coeff[3] * m +
        coeff[4] * wl +
        coeff[5] * m**2 * wl +
        coeff[6] * m * wl**2 +
        coeff[7] * m**3 * wl**3

    This version of the function to evaluate the dispersion relation
    assumes that the grating is first order, i.e. m = 1.  The dispersion
    coefficients give one-indexed pixel coordinates (reference pixels),
    but this function converts to zero-indexed pixels.

    @param coeff:  a list of eight elements containing the
        dispersion coefficients as read from a STIS _dsp.fits table
    @type coeff:  either a list of float or a numarray of Float64

    @param wl:  a single wavelength or an array (numarray) of wavelengths,
        in Angstroms
    @type wl:  either a float or a numarray of Float64

    @return:  the pixel number (or array of pixel numbers) corresponding
        to the input wavelength(s); note that these are zero indexed
    @rtype:  either a float or a numarray of Float64
    """

    c = [0., 0., 0., 0.]
    c[0] = coeff[0] + coeff[3]
    c[1] = coeff[1] + coeff[4] + coeff[5]
    c[2] = coeff[2] + coeff[6]
    c[3] = coeff[7]

    # x = c[0] + c[1] * wl + c[2] * wl**2 + c[3] * wl**3
    x = c[3] * wl
    x = (c[2] + x) * wl
    x = (c[1] + x) * wl
    x = c[0] + x

    return x - 1.       # zero indexed
