

def newton(x, coeff, cenwave, niter=4):
    """Return the wavelength corresponding to pixel x.

    The dispersion solution is evaluated iteratively, and the slope
    (dispersion) for Newton's method is determined numerically, using a
    difference in wavelength of one Angstrom.  Note that the evalDisp
    in this file assumes that the grating is first order.

    Parameters
    -----------
    x : float or ndarray
        a single pixel number or an array of pixel numbers
    coeff : array_like object
        a list of eight elements containing the
        dispersion coefficients as read from a STIS _dsp.fits table
    cenwave : int or float
        central wavelength, in Angstroms
    niter : int
        number of iterations

    Returns
    -------
    wavelength : float or ndarray
        a single wavelength or an array (numarray) of wavelengths,
        in Angstroms
    """

    wl = cenwave
    x0 = evalDisp(coeff, wl)

    delta_wl = 1.                       # one Angstrom
    for i in range(niter):
        x1 = evalDisp(coeff, wl+delta_wl)
        dispersion = delta_wl / (x1 - x0)
        wl += dispersion * (x - x0)
        x0 = evalDisp(coeff, wl)

    return wl


def evalDisp(coeff, wl):
    """Return the pixel corresponding to wavelength wl.

    Notes
    -----
    The expression in the calstis code is::

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

    Parameters
    -----------
    coeff : array_like object
        a list of eight elements containing the
        dispersion coefficients as read from a STIS _dsp.fits table

    wl : float or ndarray
        a single wavelength or an array (numarray) of wavelengths,
        in Angstroms

    Returns
    --------
    pix_number : float or ndarray
        the pixel number (or array of pixel numbers) corresponding
        to the input wavelength(s); note that these are zero indexed
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
