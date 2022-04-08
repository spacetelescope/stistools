from astropy.io import fits
import numpy as np

__doc__ = """
The purpose of this ctestis task is to correct signal levels of point-like
sources in photometry tables measured from STIS CCD images for charge loss
due to imperfect Charge Transfer Efficiency (CTE). The algorithm used to
correct for CTE-induced signal loss is based on the equations published in 
Goudfrooij, Bohlin, Maiz-Apellaniz, & Kimble, PASP, October 2006 edition
(astro-ph/0608349). The values of CTE loss derived using this algorithm
should be accurate to about 3% RMS (tested for data taken between March
1997 and August 2004). No significant differences in CTE loss were found
for different aperture sizes, although this has been verified only for a
limited range of aperture sizes (2, 3, and 5 pixel radii). The algorithm
was derived from measurements of point sources in a relatively sparse
field (the outskirts of a Galactic globular cluster), as detailed in the
PASP paper mentioned above.

The function also computes the shift in the Y centroid of point sources
due to distortions in the steller PDF cause by CTE trails. the algorithm
is taken from the Equation 9 of Goodfrooij et al. (2006). Note, however,
that the equation has been multiplied by -1, so that the resulting
correction may be ADDED to measured Y centriod of each star.

The code takes inputs of net counts for a source (background subtracted), a
sky-background estimate, and the source Y-position on the detector (since
CTI effects worsen furthest from the readout). The single-pixel
sky-background estimate should be measured from individual cosmic-ray (CR)
split, bias- and dark-subtracted, and flat-fielded images (flt.fits) that
have not had any sky subtracted. This can determined with random sampling
and/or iterative sigma-clipping of sky pixels (e.g., Goudfrooij et al.
2006). The net counts measured from the science images (summed, sky-
subtracted exposures) should then be scaled to the exposure time of the CR
split FLT image (e.g., if CRSPLIT=5, the net counts are divided by five).
Note that not all extensions composing an FLT file necessarily have equal
exposure times, so a fractional scaling of the CR split to total exposure
time of the CR-combined science image (e.g. CRJ) science image should be
used to scale the counts. The magnitude correction (dmagc) measured using
single-CR split parameters can be added to the magnitude derived from the
total exposure time science image with no further scaling.

If working with CRSPLIT scaled sky and net counts values, the filename 
(stisimage) should not be provided to avoid pulling incorrect information 
from the image headers. The following parameters therefore should be set 
manually: nread=1 (indicating it is just one CRSPLIT exposures), gain, 
mjd (start date), and ybin (BINAXIS2). If image is an sx2.fits file, set 
sx2=True.


Examples
--------

:func:`ctestis` with ycol set to 182, net set to 5,000 and sky set to 150.

>>> from stistools.ctestis import ctestis
>>> fluxc, dmagc, dyc = ctestis(182., 5000., 150., stisimage='o4qp9g010_crj.fits')
mjd: 50893.30
nread: 2
ybin: 1
gain:  1.0
amp: D
tt0: -2.3865942
lcts: -0.67595399
bck: 75.0
lbck: 2.317577
cti: 1.7314006e-05
fluxc: 2536.7133
dmagc: -0.015828427
cti10000: 0.17314006
dy512: 0.0043051192
dyc: 0.007079903
net: 5000.0
sky: 150.0
ycol: 182.0
fluxc: 2536.7133
dmagc: -0.015828427
dyc: 0.007079903

"""

__taskname__ = "ctestis"
__version__ = "1.0"
__vdate__ = "25-January-2019"
__author__ = "Python version (2018): Sara Ogaz, " \
             "IDL version (2015) : Sean Lockwood (edits from John Biretta), " \
             "CL version (2006) : P. Goudfrooij (edits from V. Dixon)"


def ctestis(ycol, net, sky, stisimage=None, mjd=None, nread=None,
                  ybin=None, gain=None, amp='D', sx2=False):
    """
    Calculate the STIS empirical correction to magnitude and astrometric shift,
    given photometry results.

    Parameters
    ----------
    ycol : arr
        Y-column # of object
    net : arr
        Net photometric counts (background subtracted) measured from science image
        and scaled to cosmic-ray split exposure time (from which sky is measured)
    sky : arr
        Single-pixel sky-background estimate measured from individual cosmic-ray
        split, bias- and dark-subtracted, flat-fielded images (flt.fits) with no 
        sky subtraction
    stisimage : str, optional
        The name of the SX2 file from which to pull the header keywords.
    mjd : float, optional
        Modified julian date corresponding to the start time of the 1st
        exposure, corresponds to the TEXPSTRT keyword. If stisimage file is
        defined TEXPSTRT keyword will overwrite any provided mjd value.
    nread : int, optional
        Number of image sets combined during CR rejection, corresponds to the
        NCOMBINE keyword. If stisimage file is defined NCOMBINE keyword will
        overwrite any provided nread value.
    ybin : int , optional
        Axis2 data bin size in unbinned detector pixels, corresponds to the
        BINAXIS2 keyword. If stisimage file is defined BINAXIS2 keyword will
        overwrite any provided ybin value.
    gain : float, optional
        The image gain, corresponds to the CCDGAIN keyword.  If stisimage file
        is defined the CCDGAIN keyword will overwrite any provided gain value.
        If the gain is 4.0, it will be upated to 4.08.
    amp : str, optional
        The amplifier used for the observation (default 'D').  Ignored if
        stisimage is provided.
    sx2 : bool, optional
        Force the procedure to remove the top/bottom 38 rows. This is
        automatically done if the file in stisname contains '_sx2'. Default
        values is False

    Returns
    -------
    fluxc : arr
        The empirically-corrected flux (counts)
    dmagc : arr
        The empirical photometric correction (delta mag)
    dyc : arr
        The empirical astrometric correction (delta pixels)
    """

    # convert any iterable input to a numpy array
    ycol = np.asarray(ycol)
    net = np.asarray(net)
    sky = np.asarray(sky)


    if stisimage is not None:
        hdulist = fits.open(stisimage)

        mjd = hdulist[0].header['TEXPSTRT']
        nread = hdulist[1].header['NCOMBINE']
        ybin = hdulist[0].header['BINAXIS2']
        gain = hdulist[0].header['CCDGAIN']
        amp = hdulist[0].header['CCDAMP'].strip().upper()

        # not sure if I really need this rootname splitter or not

    else:
        opt_inputs = {'mjd': mjd, 'nread': nread, 'ybin': ybin, 'gain': gain}
        for key, eleme in opt_inputs.items():
            if eleme is None:
                raise ValueError("If no image filename is specified, mjd, "
                                 "nread, ybin, and gain must be specified. {} "
                                 "is missing".format(key))

        amp = amp.upper()

        # check that amp makes sense
        if amp not in ['A', 'B', 'C', 'D']:
            raise ValueError("Amplifier must be 'A', 'B', 'C', or 'D'. {} is "
                             "note a recognized value".format(amp))

        print("MJD   = ", mjd)
        print("NREAD = ", nread)
        print("YBIN  = ", ybin)
        print("GAIN  = ", gain, " --> 4.08" if gain == 4.0 else "")

    if gain == 4.0:
        gain = 4.08

    # cte equation constants from table 8
    a = 0.000133
    b = 0.54
    c = 0.205
    d = 0.05
    e = 0.82
    f = 3.6
    g = 0.21

    # equation 8 inputs
    tt0 = (mjd - 51765) / 365.25
    cts = np.maximum((net * gain / nread), 1)
    bck = np.maximum((sky * gain / nread), 0)
    lcts = np.log(cts) - 8.5
    lbck = np.log(np.sqrt(bck * bck + 1)) - 2

    # equation 8 for image cti
    cti1 = a * np.exp(-b * lcts) * (c * tt0 + 1)
    cti2 = d * np.exp(-e * lbck) + (1 - d) * np.exp(-f * ((bck / cts) ** g))
    cti = cti1 * cti2

    # equation 9 for the shift at the central row
    cti10000 = cti * 10000
    dy512 = 0.025 * cti10000 - (0.00078 * cti10000 * cti10000)

    # prep ycol
    if sx2 or ((stisimage is not None) and ("_sx2" in stisimage)):
        remove_buffer = 38
    else:
        remove_buffer = 0

    ycol -= remove_buffer  # remove the bottom 38 pixels from _SX2 files

    # Double check direction for other amplifiers
    if amp in ['B', 'D']:
        ycol_dir = ycol
    else:
        ycol_dir = 1024. / ybin - ycol

    # Outputs
    # equation 11 for corrected counts
    fluxc = cts / ((1 - cti) ** (1024. - ycol_dir * ybin))
    # mag correction
    dmagc = 2.5 * np.log10(cts / fluxc)
    # scale central row shift for other rows
    dyc = dy512 * ((1024. - ycol_dir) / 512)

    # print results

    print("")
    print("mjd: {:8.2f}\n"
          "nread: {:6}\n"
          "ybin: {:6}\n"
          "gain: {:4.2}\n"
          "amp: {}\n".
          format(mjd, str(nread), str(ybin), float(gain), amp))

    if type(lcts) is np.ndarray:
        print("tt0: {:.8}\n"
              "lcts: {}\n"
              "bck: {}\n"
              "lbck: {}\n"
              "cti: {}\n"
              "fluxc: {}\n"
              "dmagc: {}\n"
              "cti10000: {}\n"
              "dy512: {}\n"
              "dyc: {}\n".
              format(tt0, lcts, bck, lbck, cti, fluxc, dmagc, cti10000,
                    dy512, dyc))
        print('net: {}\n'
              'sky: {}\n'
              'ycol: {}\n'
              'fluxc: {}\n'
              'dmagc: {}\n'
              'dyc: {}\n'.format(net, sky, ycol, fluxc, dmagc, dyc))

    else:
        print("tt0: {:.8}\n"
              "lcts: {:.8}\n"
              "bck: {:.8}\n"
              "lbck: {:.8}\n"
              "cti: {:.8}\n"
              "fluxc: {:.8}\n"
              "dmagc: {:.8}\n"
              "cti10000: {:.8}\n"
              "dy512: {:.8}\n"
              "dyc: {:.8}\n".
              format(tt0, lcts, bck, lbck, cti, fluxc, dmagc, cti10000,
                     dy512, dyc))

        print('net: {:.8}\n'
              'sky: {:.8}\n'
              'ycol: {:.8}\n'
              'fluxc: {:.8}\n'
              'dmagc: {:.8}\n'
              'dyc: {:.8}\n'.format(net, sky, ycol, fluxc, dmagc, dyc))

    return fluxc, dmagc, dyc
