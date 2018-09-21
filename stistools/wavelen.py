
import numpy as N

from . import evaldisp
from . import gettable
from . import radialvel
from . import r_util

DEG_RAD = N.pi / 180.                   # degrees to radians
SPEED_OF_LIGHT = 299792.458             # km / s


def compute_wavelengths(shape, phdr, hdr, helcorr):
    """Compute a 2-D array of wavelengths, one value for each image pixel.

    Parameters
    ----------
    shape : tuple of two ints
        the number of rows and columns in the output image
    phdr :  fits Header object
        primary header
    hdr : fits Header object
        extension header
    helcorr : string
        "PERFORM" if heliocentric correction should be done

    Returns
    -------
    wavelengths : ndarray of float64
        an array of wavelengths, of the same shape (nrows, ncols) as
        the output image

    """

    REF_ANGLE = 0.315                           # degrees
    PSEUDOAPERTURES = ["E1", "E2", "D1"]        # currently defined pseudo-aps

    (nrows, ncols) = shape

    opt_elem = phdr.get("opt_elem", default="keyword_missing")
    cenwave = phdr.get("cenwave", default=0)
    aperture = phdr.get("aperture", default="keyword_missing")
    propaper = phdr.get("propaper", default="keyword_missing")
    sclamp = phdr.get("sclamp", default="NONE")
    disptab = phdr.get("disptab", default="keyword_missing")
    apdestab = phdr.get("apdestab", default="keyword_missing")
    inangtab = phdr.get("inangtab", default="keyword_missing")
    ra_targ = phdr.get("ra_targ", default="keyword_missing")
    dec_targ = phdr.get("dec_targ", default="keyword_missing")

    expstart = hdr.get("expstart", default=0.)
    expend = hdr.get("expend", default=0.)
    crpix2 = hdr.get("crpix2", default=0.)
    ltm = hdr.get("ltm1_1", default=1.)
    ltv1 = hdr.get("ltv1", default=0.)
    ltv2 = hdr.get("ltv2", default=0.)
    shifta1 = hdr.get("shifta1", default=0.)
    shifta2 = hdr.get("shifta2", default=0.)

    disptab = r_util.expandFileName(disptab)
    apdestab = r_util.expandFileName(apdestab)
    inangtab = r_util.expandFileName(inangtab)

    # Modify ltv and crpix2 for zero-indexed pixels.
    ltv1 += (ltm - 1.)
    crpix2 -= 1.
    binaxis1 = round(1. / ltm)         # should be 1, 2, 4 or 8

    # These offsets have not been converted from one-indexing to
    # zero-indexing, but that's OK because the input image must not
    # be binned in the cross-dispersion (axis 2) direction.
    offset = shifta2 - ltv2

    if sclamp != "NONE":
        nchar = len(propaper)
        ending = propaper[nchar-2:nchar]
        if ending in PSEUDOAPERTURES:
            aperture = aperture + ending

    if helcorr == "PERFORM":
        v_helio = radialvel.radialVel(ra_targ, dec_targ,
                                      (expstart + expend) / 2.)
        hfactor = (1. - v_helio / SPEED_OF_LIGHT)
    else:
        hfactor = 1.

    # Get the dispersion relation.  This will be in the form of
    # a set of coefficients at several positions along the slit, i.e.
    # disp_coeff[i] is the set of coefficients at position a2center[i].
    filter = {"opt_elem": opt_elem, "cenwave": cenwave}
    disp_info = gettable.getTable(disptab, filter,
                                  sortcol="a2center", at_least_one=True)
    ref_aper = disp_info.field("ref_aper")[0]  # name of reference aperture
    a2center = disp_info.field("a2center") - 1.        # zero indexing
    ncoeff = disp_info.field("ncoeff")[0]      # same for all rows
    disp_coeff = disp_info.field("coeff")
    delta_offset1 = get_delta_offset1(apdestab, aperture, ref_aper)

    apdes_info = gettable.getTable(apdestab, {"aperture": aperture},
                                   exactly_one=True)
    # Check whether ANGLE is a column in this table.
    names = []
    for name in apdes_info.names:
        names.append(name.lower())
    if "angle" in names:
        angle = apdes_info.field("angle")[0]
    else:
        print("Warning:  Column ANGLE not found in", apdestab)
        angle = REF_ANGLE
    del names

    delta_tan = N.tan(angle * DEG_RAD) - N.tan(REF_ANGLE * DEG_RAD)

    # Note:  this assumes a first-order spectrum, but at the time of
    # writing there's actually no distinction in any of the iac tables.
    filter = {"opt_elem": opt_elem, "cenwave": cenwave, "sporder": 1}
    inang_info = gettable.getTable(inangtab, filter, exactly_one=True)

    wavelengths = N.zeros((nrows, ncols), dtype=N.float64)
    image_pixels = N.arange(ncols, dtype=N.float64)
    # Convert from image pixels to reference pixels (but zero indexed).
    pixels = (image_pixels - ltv1) * binaxis1
    for j in range(nrows):
        row = float(j) + offset        # account for possible subarray
        # Interpolate to get dispersion relation for current (0-indexed) row.
        coeff = r_util.interpolate(a2center, disp_coeff, row)
        # Apply corrections.
        adjust_disp(ncoeff, coeff, delta_offset1, shifta1, inang_info,
                    delta_tan, row-crpix2, binaxis1)
        # Compute wavelength from pixel numbers.
        wl = evaldisp.newton(pixels, coeff, cenwave)
        wl *= hfactor
        wavelengths[j] = wl.copy()

    return wavelengths


def get_delta_offset1(apdestab, aperture, ref_aper):
    """Get the incidence angle offset.

    Parameters
    ----------
    apdestab : string
        name of the aperture description table
    aperture : string
        aperture (slit) name
    ref_aper : string
        name of the reference aperture, the one that was
        used to calculate the dispersion relation

    Returns
    -------
    angle : float
        incidence angle offset in degrees

    """

    # Get the offset for the aperture that was used for the observation.
    apdes_info = gettable.getTable(apdestab, {"aperture": aperture},
                                   exactly_one=True)
    aperture_offset1 = apdes_info.field("offset1")[0]

    # Get the offset for the aperture that was used for creating the
    # dispersion relation.
    apdes_info = gettable.getTable(apdestab, {"aperture": ref_aper},
                                   exactly_one=True)
    ref_aper_offset1 = apdes_info.field("offset1")[0]

    return aperture_offset1 - ref_aper_offset1


def adjust_disp(ncoeff, coeff, delta_offset1, shifta1, inang_info,
                delta_tan, delta_row, binaxis1):
    """Adjust the dispersion coefficients.

    The changes to the coefficients are for the incidence angle
    correction, the offset from the SHIFTA1 keyword, and the tilt
    of the slit.  The coefficients will be modified in-place.

    Parameters
    ----------
    ncoeff : int
        number of dispersion coefficients
    coeff : ndarray of float64
        array of dispersion coefficients, modified in-place
    delta_offset1 : float
        incidence angle offset in degrees
    shifta1 : float
        MSM offset (ref. pixels) in the dispersion direction
    delta_tan : float
        difference in tangents of slit angle and ref angle
    delta_row : float
        difference between current row number and CRPIX2
    binaxis1 : float
        binning factor in dispersion direction
    inang_info : rec_array
        rows from the incidence-angle table

    """

    iac_ncoeff1 = inang_info.field("ncoeff1")[0]
    iac_coeff1 = inang_info.field("coeff1")[0]
    iac_ncoeff2 = inang_info.field("ncoeff2")[0]
    iac_coeff2 = inang_info.field("coeff2")[0]

    for i in range(iac_ncoeff1):
        coeff[i] += iac_coeff1[i] * delta_offset1

    if iac_ncoeff2 > 0:
        coeff[0] += iac_coeff2[0] * delta_offset1
    if iac_ncoeff2 > 1:
        coeff[1] += iac_coeff2[1] * delta_offset1**2

    # Correct for MSM shift.
    coeff[0] += shifta1

    # Correct for slit tilt.
    coeff[0] += (delta_tan * delta_row * binaxis1)
