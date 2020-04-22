import numpy as np
from astropy.modeling import models, fitting

# from lines 280 through 344 of mkfringeflat.cl
def find_loc(input, low_frac=0.2, high_frac=0.8, low_line_frac=0.4):
    """Find the cross-dispersion location of the target spectrum.

    Parameters
    ----------
    input: ndarray
        The input science data array for the current image set.

    low_frac: float
        Fraction of image width for the start of a slice.

    high_frac: float
        Fraction of image width for the end of a slice.

    low_line_frac: float
        Fraction of image height for the start of a slice for limiting the
        region over which to search for the maximum of the slit profile.

    Returns
    -------
    target_locn: float or None
        The location (zero based pixel coordinate) of the target in the
        cross-dispersion direction.
        None will be returned if the quadratic fit to the cross-dispersion
        profile has zero curvature.
    """

    shape = input.shape
    # first column in a slice
    fcol = int(round(low_frac * shape[1]))
    # last column in a slice
    lcol = int(round(high_frac * shape[1])) + 1

    # first line in a slice
    fline = int(round(low_line_frac * shape[0]**2 / 1024.)) + 1
    # last line in a slice
    lline = shape[1] - fline

    slit_prof = input[:, fcol:lcol].mean(axis=1, dtype=np.float64)
    subset_slit_prof = slit_prof[fline:lline]
    med = np.median(subset_slit_prof)
    sigma = subset_slit_prof.std()
    low_lim = med + 3. * sigma

    mask = np.where(subset_slit_prof >= low_lim)
    # The independent variable is the pixel index in the cross-dispersion
    # direction (the vertical direction).
    x_array = mask[0]
    # The dependent variable is the slit profile data value.
    y_array = subset_slit_prof[mask]

    # IMHO, this is not sensible, but it's what mkfringeflat.cl does.
    # tcalc ("temp_slitprof.tab", "data", \
    #        "if data < 0.0 then data * (-1.) else data", datatype="real")
    y_array = np.where(y_array < 0., -y_array, y_array)

    index = y_array.argsort()
    # The zero point for the index array is fline.
    maxrow = x_array[index[-1]] + fline

    ffrow = maxrow - 2
    flrow = maxrow + 3

    fit = fitting.LinearLSQFitter()
    x = np.arange(len(slit_prof), dtype=np.float64)
    quad = fit(models.Polynomial1D(2), x[ffrow:flrow], slit_prof[ffrow:flrow])
    if quad.c2.value == 0.:
        target_locn = None
    else:
        target_locn = -0.5 * (quad.c1.value / quad.c2.value)

    return target_locn
