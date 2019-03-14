#! /usr/bin/env python

from math import modf, sqrt
import os
import argparse

from astropy.io import fits
import numpy as np

__doc__ = """
Analyze STIS target acquisition images. :func:`tastis` will print general
information about each input target acquisition image, and will analyze both
types of STIS target acquisitions: ACQs and ACQ/PEAKs

ACQ procedure is described in "STIS Onboard CCD Target Acquisitions" in the
STIS Instrument Handbook.  The ACQ/PEAK procedure is described in "Onboard
Target-Acquisition Peakups (ACQ/PEAK)" also in the STIS Instrument Handbook.

Target positions in global and local (subarray) coordinates and the total flux
of the target in the maximum checkbox during both acquistions phases (course
and fine) are displayed.

If update=True, keywords are added to the header to make problems easier to
locate in batch processing.  Warnings are given if the spt file is not present
when :func:`tastis` is run.

Examples
--------

:func:`tastis` with the default of update=False:

>>> import stistools
>>> stistools.tastis.tastis("ocmv0lw6q_raw.fits")
===============================================================================
oc7w11viq       HST/STIS    G430L        0.3X0.05ND             ACQ/PEAK-UP
prop: 13465      visit: 11    line: 3   target: HD128621-2
obs date, time: 2014-07-24    22:05:06   exposure time:  0.10
dom GS/FGS: S7QX000330F1    sub-dom GS/FGS: S7QX000694F2
ACQ params:     bias sub: 1510                     method: RETURN-TO-BRIGHTEST
subarray (axis1,axis2):   size=(1022,32)          corner=(25,500)
-------------------------------------------------------------------------------
Scan type: LINEARAXIS2                  Step size (mas): 250
 [210 753   0]
                    axis1 axis2     axis1  axis2         V2   V3
                      (pixels)        (arcsec)           (arcsec)
Estimated slew:      0.0  -0.1     0.000 -0.005     -0.004  0.004
Flux in post-slew confirmation image (751752) - Pedestal (748587) = 3165 DN
-------------------------------------------------------------------------------
The flux in the confirmation image is 320% greater than the maximum flux
in the ACQ/PEAK scan.  An excess greater than 100% indicates
problems in the ACQ/PEAK.
The flux in the confirmation image is 16% of the recommended minimum
of 20000 DN for a dispersed-light ACQ/PEAK.  The signal-to-noise in
the AC
===============================================================================

:func:`tastis` with update=True:

>>> import stistools
>>> stistools.tastis.tastis("ocmv0lw6q_raw.fits", update=True)
===============================================================================
ocmv0lw6q       HST/STIS    MIRVIS      F25ND3             ACQ/POINT
prop: 13760      visit: 0L    line: 1   target: CD-59D3300
obs date, time: 2016-09-29    23:43:50   exposure time:  1.10
dom GS/FGS: S4B0000993F2    sub-dom GS/FGS: S4B0000953F1
ACQ params:     bias sub: 1510   checkbox: 3      method: FLUX CENTROID
subarray (axis1,axis2):   size=(100,100)          corner=(487,466)
-------------------------------------------------------------------------------
Coarse locate phase:           Target flux in max checkbox (DN): 1560
                       global          local
                    axis1 axis2     axis1 axis2
Target location:    534.2  507.0    48.2  42.0
                    axis1 axis2     axis1  axis2         V2      V3
                      (pixels)        (arcsec)            (arcsec)
Estimated slew:     -1.5  -9.0    -0.079 -0.457       -0.379  0.268
-------------------------------------------------------------------------------
Fine locate phase:            Target flux in max checkbox (DN): 1559
                       global            local
                    axis1 axis2     axis1 axis2
Target location:    534.2  516.8    48.2  51.8
Ref ap location:    537.5  517.0    19.5  17.0
                    axis1 axis2     axis1  axis2         V2      V3
                      (pixels)        (arcsec)           (arcsec)
Estimated slew:     -2.1  -0.2     -0.104 -0.010      -0.081 -0.067
-------------------------------------------------------------------------------
Total est. slew:    -3.6  -9.2    -0.183 -0.467        -0.460  0.201
-------------------------------------------------------------------------------
Your ACQ appears to have succeeded, as the fluxes in the coarse
and fine stages agree within 25% and the fine slews were less than
4 pixels as expected
===============================================================================

"""

__taskname__ = "tastis"
__version__ = "1.0"
__vdate__ = "14-June-2018"
__author__ = "Python: Sara Ogaz.  C code: R. Katsanis, R. Downes, " \
             "Phil Hodge. Original IDL code: George Hartig"

# These are possible values for the bit mask badacq.
BAD_ACQ = 1	            # any problem
BAD_SLEW = 2	        # ACQ only
BAD_LAMP_LOW = 4	    # ACQ only
BAD_RATIO_HIGH = 8	    # ACQ or ACQ/PEAK
BAD_RATIO_LOW = 16	    # ACQ or ACQ/PEAK
BAD_SATURATED = 32	    # ACQ or ACQ/PEAK
BAD_FLUX = 64	        # ACQ/PEAK only
BAD_END = 128	        # ACQ/PEAK only
BAD_TDF = 256	        # take-data flag was down

FATAL_ERROR = 2
LOW_FLUX_CUTOFF = 0.8
HIGH_FLUX_CUTOFF = 2.0
MAX_GOODMAX = 32000.	 # higher would be saturated */
MIN_GOODMAX = 1900.	     # lower implies lamp was not on */
MIN_IMAGING_FLUX = 1250.
MIN_SPECTROSCOPIC_FLUX = 20000.

PLATESCALE = 0.0508
COSTHETA = 0.70711
SINTHETA = 0.70711


def tastis(raw_filename, update=False):
    """
    Analyze STIS target acquisition images.

    Parameters
    ----------
    raw_filename: str
        Name of the input raw file. For some raw files you will need a copy of
        the spt file in the same directory.

    update: bool
        If True, keywords associated with tastis checks will be updated.
        Default values is False.
    """

    # filename can take wildcards, worry about this later
    # checks for matching files exits if it can't find any
    # loops over files

    # build string for raw and spt file
    spt_filename = raw_filename.replace("raw", "spt")
    spt_exists = os.path.exists(spt_filename)

    # open file
    head = fits.getheader(raw_filename)
    obsmode = head['OBSMODE']

    if obsmode == "ACQ" or obsmode == "ACQ/PEAK":

        # Check for spt file if obsmode is ACQ/PEAK
        if not spt_exists and obsmode == "ACQ/PEAK":
            FileNotFoundError("Can't find {} (required for ACQ/PEAK), "
                              "exiting".format(spt_filename))

    keywords = _read_keywords(raw_filename, spt_filename, spt_exists)

    _calculate_slews(keywords)

    # Read dominant & subdominant FGS from keywords DGESTAR & SGESTAR
    # in primary header of spt file.
    if spt_exists:
        spt_head = fits.getheader(spt_filename)
        keywords['domfgs'] = spt_head['dgestar']
        keywords['subfgs'] = spt_head['sgestar']
    else:
        keywords['domfgs'] = ""
        keywords['subfgs'] = ""

    badacq = _print_output(keywords, spt_exists)

    # Update keywords in the input primary header to indicate
    # which tests succeeded and which failed.
    if update:
        with fits.open(raw_filename, mode='update') as raw_hdulist:
            if badacq:
                raw_hdulist[0].header['acqstat'] = "FAILED"
            else:
                raw_hdulist[0].header['acqstat'] = "OK"

            if keywords['obsmode'] == "ACQ":
                # Bad ratio
                if badacq & BAD_RATIO_HIGH:
                    raw_hdulist[0].header['acq_rat'] = "HIRATIO"
                elif badacq & BAD_RATIO_LOW:
                    raw_hdulist[0].header['acq_rat'] = "LORATIO"
                else:
                    raw_hdulist[0].header['acq_rat'] = "OKRATIO"

                # Bad slew
                if badacq & BAD_SLEW:
                    raw_hdulist[0].header['acq_slew'] = "BIGSLEW"
                else:
                    raw_hdulist[0].header['acq_slew'] = "OK_SLEW"

                # Saturation
                if badacq & BAD_SATURATED:
                    raw_hdulist[0].header['acq_sat'] = "SAT"
                else:
                    raw_hdulist[0].header['acq_sat'] = "UNSAT"

                # Bad Lamp
                if badacq & BAD_LAMP_LOW:
                    raw_hdulist[0].header['acq_lamp'] = "LO_LAMP"
                else:
                    raw_hdulist[0].header['acq_lamp'] = "OK_LAMP"

            # for ACQ/PEAK
            else:
                # Bad ratio
                if badacq & BAD_RATIO_HIGH:
                    raw_hdulist[0].header['acqp_rat'] = "HIRATIO"
                elif badacq & BAD_RATIO_LOW:
                    raw_hdulist[0].header['acqp_rat'] = "LORATIO"
                else:
                    raw_hdulist[0].header['acqp_rat'] = "OKRATIO"

                # Bad flux
                if badacq & BAD_FLUX:
                    raw_hdulist[0].header['acqp_flx'] = "LO_FLUX"
                else:
                    raw_hdulist[0].header['acqp_flx'] = "OK_FLUX"

                # Saturation
                if badacq & BAD_SATURATED:
                    raw_hdulist[0].header['acqp_sat'] = "SAT"
                else:
                    raw_hdulist[0].header['acqp_sat'] = "UNSAT"

                # Bad end
                if badacq & BAD_END:
                    raw_hdulist[0].header['acqp_end'] = "HI_END"
                else:
                    raw_hdulist[0].header['acqp_end'] = "OK_END"

            # Bad tdf
            if badacq & BAD_TDF:
                raw_hdulist[0].header['dataflag'] = "TDFDown"
            elif not spt_exists:
                raw_hdulist[0].header['dataflag'] = "UNKNOWN"
            else:
                raw_hdulist[0].header['dataflag'] = "TDF_Up"


def _read_keywords(raw_filename, spt_filename, spt_exists):
    """
    Read in raw and spt FITS file header keywords used in :func:`tastis`, and
    store the results in a dictionary returned by the function.

    Parameters
    ----------
    raw_filename: str
        Name of the FITS input raw file.

    spt_filename: str
        Name of the FITS spt file associated with the raw file.

    spt_exists: bool
        If True, the `spt_filename` exists in the same directory as the raw
        file.

    Returns
    -------
    keywords: dict
        dictionary containing all keywords and other data needed by
        :func:`tastis`.
    """

    keywords = {}

    with fits.open(raw_filename) as raw_hdulist:

        # Read universal primary header keywords
        prim_header_keywords = ['rootname', 'obsmode', 'obstype', 'proposid',
                                'sizaxis1', 'sizaxis2', 'texptime', 'biaslev',
                                'targname', 'tdateobs', 'ttimeobs', 'linenum',
                                'centera1', 'centera2']

        for key in prim_header_keywords:
            keywords[key] = raw_hdulist[0].header[key]

        keywords['optelem'] = raw_hdulist[0].header['opt_elem']

        # For the aperture name, use PROPAPER if its value ends in "E1"
        # or "D1"; otherwise, use APERTURE.
        propaper = raw_hdulist[0].header['propaper']
        if propaper[-2:] in ['E1', 'D1']:
            keywords['aperture'] = propaper
        else:
            keywords['aperture'] = raw_hdulist[0].header['aperture']

        # grab from 1st ext (but not into dict) ngoodpix, goodmean
        ngoodpix = raw_hdulist[1].header['ngoodpix']
        goodmean = raw_hdulist[1].header['goodmean']

        # Obsmode dependent header pulls
        if keywords['obsmode'] == 'ACQ':
            # 0th header
            keywords['acqtype'] = raw_hdulist[0].header['acqtype']
            keywords['box_step'] = raw_hdulist[0].header['checkbox']

            # 1th header
            keywords['counts1'] = raw_hdulist[1].header['maxchcnt']
            keywords['targax1'] = raw_hdulist[1].header['targa1']
            keywords['targay1'] = raw_hdulist[1].header['targa2']

            # 4th header
            keywords['counts2'] = raw_hdulist[4].header['maxchcnt']
            keywords['goodmax2'] = raw_hdulist[4].header['goodmax']
            keywords['targax4'] = raw_hdulist[4].header['targa1']
            keywords['targay4'] = raw_hdulist[4].header['targa2']

            # 7th header
            keywords['goodmax3'] = raw_hdulist[7].header['goodmax']
            keywords['apera1'] = raw_hdulist[7].header['apera1']
            keywords['apera2'] = raw_hdulist[7].header['apera2']
            keywords['aperlka1'] = raw_hdulist[7].header['aperlka1']
            keywords['aperlka2'] = raw_hdulist[7].header['aperlka2']

            # set pedestal and goodmax1 to 0, not used for "ACQ"
            keywords['pedestal'] = 0
            keywords['goodmax1'] = 0

            if keywords['acqtype'] == "POINT":
                keywords['search'] = "FLUX CENTROID"
            else:
                keywords['search'] = raw_hdulist[0].header['centmeth']

        # Read keywords from the ACQ/PEAK primary header and from the
        # spt extension header.
        if keywords['obsmode'] == "ACQ/PEAK":
            keywords['peakcent'] = raw_hdulist[0].header['peakcent']
            keywords['search'] = raw_hdulist[0].header['pksearch']
            keywords['box_step'] = raw_hdulist[0].header['numsteps']
            keywords['peakstep'] = raw_hdulist[0].header['peakstep']
            keywords['pedestal'] = raw_hdulist[0].header['pedestal']

            keywords['goodmax1'] = raw_hdulist[1].header['goodmax']

            # From spt file 1st header
            spt_head = fits.getheader(spt_filename, ext=1)
            keywords['otaslwa1'] = spt_head['otaslwa1']
            keywords['otaslwa2'] = spt_head['otaslwa2']

            # Calculate post-slew flux, get pedestal (from raw file) & dwell
            # fluxes. For ACQ/PEAKs only.
            keywords['counts2'] = 0.  # not used for ACQ/PEAK
            keywords['goodmax2'] = 0.  # not used for ACQ/PEAK
            keywords['goodmax3'] = 0.  # not used for ACQ/PEAK
            keywords['counts1'] = ngoodpix * goodmean

            # Read dwell fluxes from 4th extension of raw file.
            # I think this whole section (lines 564-573) are just pulling
            # the data array from the 4th image extension (raw file)
            keywords['naxis1'] = raw_hdulist[4].data.shape[1]
            keywords['naxis2'] = raw_hdulist[4].data.shape[0]
            keywords['dwell'] = raw_hdulist[4].data

        raw_hdulist.close()

    # check for spt file
    if spt_exists:
        spt_head = fits.getheader(spt_filename, ext=1)
        keywords['ocstdfx'] = spt_head['ocstdfx']
    else:
        keywords['ocstdfx'] = "unknown"

    # Extract visit & expnum from linenum before period goes into keyword dict
    # as visit. if this is the end of the string, fill expnum in the keyword
    # dict with 0, otherwise fill expnum with rest of linenum value after
    # period converted to float
    split_linenum = keywords['linenum'].split(".")
    keywords['visit'] = split_linenum[0]
    if len(split_linenum) == 1:
        keywords['expnum'] = 0
    else:
        keywords['expnum'] = float(split_linenum[1])

    # Calculate corner from 'centera' & sizaxis'.
    keywords['corner1'] = keywords['centera1'] - keywords['sizaxis1']/2
    keywords['corner2'] = keywords['centera2'] - keywords['sizaxis2']/2

    # Calculate coarse, fine local axis & reference aperture locations.
    # For ACQs only.
    if keywords['obsmode'] == 'ACQ':
        keywords['coarse1'] = keywords['targax1'] - (keywords['corner1'] - 1) + 1
        keywords['coarse2'] = keywords['targay1'] - (keywords['corner2'] - 1) + 1
        keywords['fine1'] = keywords['targax4'] - (keywords['corner1'] - 1) + 1
        keywords['fine2'] = keywords['targay4'] - (keywords['corner2'] - 1) + 1
        keywords['refaper1'] = keywords['apera1'] - (keywords['corner1'] + 31) +1
        keywords['refaper2'] = keywords['apera2'] - (keywords['corner2'] + 34) +1
        if keywords['box_step'] > 3:
            offset = (keywords['box_step'] + 1)/2
            keywords['refaper1'] -= offset
            keywords['refaper2'] -= offset

    return keywords


def _calculate_slews(keywords):
    """
    Calculate slew information used by :func:`tastis` using input data
    dictionary. THIS FUNCTION EDITS THE DICTIONARY IN PLACE.

    Parameters
    ----------
    keywords: dict
        dictionary containing all keywords and other data needed for slew
        calculation. DICTIONARY IS EDITED IN PLACE.
    """

    # Slew calculation for ACQs.
    if keywords['obsmode'] == 'ACQ':

        # Define all possible apertures for ACQs.
        aperture_acq_dict = {"F25NDQ1": -1.24840,
                             "F25NDQ2": -1.24840,
                             "F25NDQ3": -1.24840,
                             "F25NDQ4": -1.24840,
                             "F28X50LP": -1.26850,
                             "F28X50OIII": -1.26850,
                             "F28X50OII": -1.31570,
                             "F25ND3": -1.24840,
                             "F25ND5": -1.24840}

        if keywords['aperture'] in aperture_acq_dict:
            offset = aperture_acq_dict[keywords['aperture']]
        else:
            offset = 0.0

        # Slews in pixels.
        keywords['a1coarse_pix'] = keywords['targax1'] - offset - keywords['aperlka1'] + 1
        keywords['a2coarse_pix'] = keywords['targay1'] - keywords['aperlka2'] + 1
        keywords['a1fine_pix'] = keywords['targax4'] - offset - keywords['apera1']
        keywords['a2fine_pix'] = keywords['targay4'] - keywords['apera2']
        keywords['a1total_pix'] = keywords['a1coarse_pix'] + keywords['a1fine_pix']
        keywords['a2total_pix'] = keywords['a2coarse_pix'] + keywords['a2fine_pix']

        # Slews in arcseconds.
        keywords['a1coarse_arc'] = _arcseconds(keywords['a1coarse_pix'])
        keywords['a2coarse_arc'] = _arcseconds(keywords['a2coarse_pix'])
        keywords['a1fine_arc'] = _arcseconds(keywords['a1fine_pix'])
        keywords['a2fine_arc'] = _arcseconds(keywords['a2fine_pix'])
        keywords['a1total_arc'] = _arcseconds(keywords['a1total_pix'])
        keywords['a2total_arc'] = _arcseconds(keywords['a2total_pix'])

        keywords['V2coarse'] = _v2coord(keywords['a2coarse_arc'],
                                        keywords['a1coarse_arc'])
        keywords['V3coarse'] = _v3coord(keywords['a1coarse_arc'],
                                        keywords['a2coarse_arc'])
        keywords['V2fine'] = _v2coord(keywords['a2fine_arc'],
                                      keywords['a1fine_arc'])
        keywords['V3fine'] = _v3coord(keywords['a1fine_arc'],
                                      keywords['a2fine_arc'])
        keywords['V2total'] = _v2coord(keywords['a2total_arc'],
                                       keywords['a1total_arc'])
        keywords['V3total'] = _v3coord(keywords['a1total_arc'],
                                       keywords['a2total_arc'])

    else:
        # Slew calculations for ACQ/PEAKs.
        if keywords['search'] == "LINEARAXIS2":
            finalx = int(keywords['box_step']/2) * \
                keywords['peakstep']/(PLATESCALE*1000.0)
            finaly = 0.0

        elif keywords['search'] == "LINEARAXIS1":
            finalx = 0.0
            finaly = int(keywords['box_step']/2) * \
                keywords['peakstep']/(PLATESCALE*1000.0)

        elif keywords['search'] == 'SPIRAL':
            x, finaly = modf(sqrt(keywords['box_step']) / 2)
            finaly = -1 * finaly * keywords['peakstep'] / (PLATESCALE*1000.0)

            x, finalx = modf(sqrt(keywords['box_step']) / 2)
            finalx = -1 * finalx * keywords['peakstep'] / (PLATESCALE*1000.0)

        # Final slews in pixels.
        keywords['a1total_pix'] = keywords['otaslwa1']/10.0 + finaly
        keywords['a2total_pix'] = keywords['otaslwa2']/10.0 + finalx

        if keywords['search'] == "SPIRAL":
            if abs(keywords['a2total_pix']) < 0.05:
                keywords['a2total_pix'] = 0.0

            if abs(keywords['a1total_pix']) < 0.05:
                keywords['a1total_pix'] = 0.0

        # Rounding up the pixel values up to the decimal place.
        keywords['a1total_pix'] = _ndec(keywords['a1total_pix'])
        keywords['a2total_pix'] = _ndec(keywords['a2total_pix'])

        # Slews in arcseconds.
        keywords['a1total_arc'] = _arcseconds(keywords['a1total_pix'])
        keywords['a2total_arc'] = _arcseconds(keywords['a2total_pix'])
        keywords['V2total'] = _v2coord(keywords['a2total_arc'],
                                       keywords['a1total_arc'])
        keywords['V3total'] = _v3coord(keywords['a1total_arc'],
                                       keywords['a2total_arc'])


def _print_output(keywords, spt_exists):
    """
    Print analysis output to stdout for :func:`tastis` report.

    Parameters
    ----------
    keywords: dict
        dictionary containing all keywords and other data needed by
        :func:`tastis`.

    spt_exists: bool
        If True, the `spt_filename` exists in the same directory as the raw
        file.

    Returns
    -------
    badacq: integer (bit flag)
        Bit flag integer containing :func:`tastis` error flags.
    """

    # Print to stdout
    print('=' * 79)

    if keywords['obsmode'] == "ACQ":
        print("{:>8}       HST/STIS    MIRVIS     {:>7}             "
              "ACQ/{}".format(keywords['rootname'], keywords['aperture'],
                              keywords['acqtype']))
    else:
        print("{:>8}       HST/STIS    {}        {:>7}             "
              "ACQ/PEAK-UP".format(keywords['rootname'], keywords['optelem'],
                                   keywords['aperture']))

    print("prop: {:4d}      visit: {}    line: {:.0f}   target: {}".format(
        keywords['proposid'], keywords['visit'], keywords['expnum'],
        keywords['targname']))

    print("obs date, time: {:>8}    {:>8}   exposure time: {:5.2f}".format(
        keywords['tdateobs'], keywords['ttimeobs'], keywords['texptime']))

    if keywords['domfgs'] != "" or keywords['subfgs'] != "":
        print("dom GS/FGS: {}    sub-dom GS/FGS: {}".
              format(keywords['domfgs'], keywords['subfgs']))

    if keywords['obsmode'] == "ACQ":
        print("ACQ params:     bias sub: {:.0f}   checkbox: {:d}      method: "
              "{}".format(keywords['biaslev'], keywords['box_step'],
                          keywords['search']))
    else:
        print("ACQ params:     bias sub: {:.0f}                     "
              "method: {}".format(keywords['biaslev'], keywords['peakcent']))

    print("subarray (axis1,axis2):   size=({:d},{:d})          "
          "corner=({:d},{:d})".format(int(keywords['sizaxis1']),
                                      int(keywords['sizaxis2']),
                                      int(keywords['corner1']),
                                      int(keywords['corner2'])))

    print('-' * 79)

    # Print rest of output according to data type: ACQ or ACQ/PEAK.
    if keywords['obsmode'] == "ACQ":
        print("Coarse locate phase:           Target flux in max checkbox "
              "(DN): {:.0f}\n".format(keywords['counts1']))
        print("                       global          local")
        print("                    axis1 axis2     axis1 axis2")
        print("Target location:    {:4.1f}  {:4.1f}    {:4.1f}  {:4.1f}\n".
              format(keywords['corner1'] + keywords['coarse1'] - 1,
                     keywords['corner2'] + keywords['coarse2'] - 1,
                     keywords['coarse1'], keywords['coarse2']))
        print("                    axis1 axis2     axis1  axis2         V2    "
              "  V3")
        print("                      (pixels)        (arcsec)            "
              "(arcsec)")
        print("Estimated slew:     {:4.1f}  {:4.1f}    {:6.3f} {:6.3f}       "
              "{:6.3f} {:6.3f}".format(keywords['a1coarse_pix'],
                                       keywords['a2coarse_pix'],
                                       keywords['a1coarse_arc'],
                                       keywords['a2coarse_arc'],
                                       keywords['V2coarse'],
                                       keywords['V3coarse']))

        # Print slews
        print('-' * 79)

        print("Fine locate phase:            Target flux in max checkbox (DN):"
              " {:.0f}\n".format(keywords['counts2']))
        print("                       global            local")
        print("                    axis1 axis2     axis1 axis2")
        print("Target location:    {:4.1f}  {:4.1f}    {:4.1f}  {:4.1f}".
              format(keywords['corner1'] + keywords['fine1'] - 1,
                     keywords['corner2'] + keywords['fine2'] - 1,
                     keywords['fine1'], keywords['fine2']))
        print("Ref ap location:    {:4.1f}  {:4.1f}    {:4.1f}  {:4.1f}\n".
              format(keywords['apera1'] + 1, keywords['apera2'] + 1,
                     keywords['refaper1'], keywords['refaper2']))
        print("                    axis1 axis2     axis1  axis2         "
              "V2      V3")
        print("                      (pixels)        (arcsec)           "
              "(arcsec)")
        print("Estimated slew:     {:4.1f}  {:4.1f}     {:6.3f} {:6.3f}      "
              "{:6.3f} {:6.3f}".format(keywords['a1fine_pix'],
                                       keywords['a2fine_pix'],
                                       keywords['a1fine_arc'],
                                       keywords['a2fine_arc'],
                                       keywords['V2fine'],
                                       keywords['V3fine']))

        # Print slews
        print('-' * 79)

        print("Total est. slew:    {:4.1f}  {:4.1f}    {:6.3f} {:6.3f}        "
              "{:6.3f} {:6.3f}". format(keywords['a1total_pix'],
                                        keywords['a2total_pix'],
                                        keywords['a1total_arc'],
                                        keywords['a2total_arc'],
                                        keywords['V2total'],
                                        keywords['V3total']))

        print('-' * 79)
        badacq = _print_warnings(keywords, spt_exists)

    else:
        print("Scan type: {}                  Step size (mas): {:.0f}".format(
            keywords['search'], keywords['peakstep']))

        if keywords['search'] == "SPIRAL":
            print("axis 1 -->,  axis 2 ^\n")

        # Print here the dwell point values
        if keywords['search'] == "LINEARAXIS2":
            # I think I might actually need to do a transpose here
            print("\n", keywords['dwell'].flatten())
        else:
            # I think I might actually need to do a transpose here (at least
            # for one of these)
            print("\n", keywords['dwell'].flatten())

        print("")

        print("                    axis1 axis2     axis1  axis2         V2   "
              "V3")
        print("                      (pixels)        (arcsec)           "
              "(arcsec)")
        print("Estimated slew:     {:4.1f}  {:4.1f}    {:6.3f} {:6.3f}     "
              "{:6.3f} {:6.3f}".format(keywords['a1total_pix'],
                                       keywords['a2total_pix'],
                                       keywords['a1total_arc'],
                                       keywords['a2total_arc'],
                                       keywords['V2total'],
                                       keywords['V3total']))
        print("Flux in post-slew confirmation image ({:.0f}) - Pedestal "
              "({:.0f}) = {:.0f} DN".
              format(keywords['counts1'], keywords['pedestal'],
                     keywords['counts1'] - keywords['pedestal']))

        print('-' * 79)
        badacq = _print_warnings(keywords, spt_exists)

    print('=' * 79)

    return badacq


def _print_warnings(keywords, spt_exists):
    """
    Print warnings output to stdout for :func:`tastis` report.

    Parameters
    ----------
    keywords: dict
        dictionary containing all keywords and other data needed by
        :func:`tastis`.

    spt_exists: bool
        If True, the `spt_filename` exists in the same directory as the raw
        file.

    Returns
    -------
    badacq: integer (bit flag)
        Bit flag integer containing :func:`tastis` error flags.
    """

    badacq = 0
    max_at_end = 0       # initial value

    if keywords['ocstdfx'] == "TDFDown":
        print("Telemetry indicates that the intended exposures may not have\n"
              "been performed.  Check the images for signal.\n")
        badacq |= BAD_TDF

    if not spt_exists:
        print("This output lacks some information because the spt.fits file\n"
              "is not present in the directory.\n")

    # ACQ warnings.
    if keywords['obsmode'] == "ACQ":
        if abs(keywords['a1fine_pix']) > 4.0 or \
                        abs(keywords['a2fine_pix']) > 4.0:
            print("The fine slew (to center the target in the reference "
                  "aperture) is larger\nthan 4 pixels.  This may indicate a "
                  "problem with your acquisition.\n")
            badacq |= BAD_SLEW

        # Ratio of flux in max checkbox in fine & coarse stages.
        ratio = keywords['counts2'] / keywords['counts1']
        if (ratio < 0.75) or (ratio > 1.25):
            print("The fluxes in the maximum checkbox in the fine and coarse "
                  "stages differ\nby more than 25%.  This may indicate a "
                  "problem with your acquisition.\n")
            if ratio < 0.75:
                badacq |= BAD_RATIO_LOW
            else:
                badacq |= BAD_RATIO_HIGH

        if keywords['goodmax2'] > MAX_GOODMAX:
            badacq |= BAD_SATURATED
            print("Saturation of pixels in the second image may have affected"
                  "\nthe final centering.\n")

        if keywords['goodmax3'] < MIN_GOODMAX:
            badacq |= BAD_LAMP_LOW
            print("The flux in the third image of the ACQ is lower than the "
                  "typical value for\n)the lamp; the image should be checked "
                  "to see if the lamp was illuminated.\n")

        if badacq == 0:
            print("Your ACQ appears to have succeeded, as the fluxes in the "
                  "coarse\nand fine stages agree within 25% and the fine "
                  "slews were less than\n4 pixels as expected\n")

    # ACQ/PEAK warnings.
    if keywords['obsmode'] == "ACQ/PEAK":
        # Calculate maximum flux in the peakup
        max_final = 0.0
        i_max = -1
        j_max = -1

        # I'm not sure if dwell would ever contain all negative values, but
        # just in case that's a possibility, to replicate original code
        # behaviour also, need to check indexing order
        if max(keywords['dwell'].flatten()) > max_final:
            max_final = max(keywords['dwell'].flatten())
            max_indexs = np.where(keywords['dwell'] ==
                                  max(keywords['dwell'].flatten()))
            i_max = max_indexs[1][0]
            j_max = max_indexs[0][0]

        # subtract pedestal
        flux = keywords['counts1'] - keywords['pedestal']
        flux_ratio = flux / max_final

        if flux_ratio < LOW_FLUX_CUTOFF:
            print("The flux in the confirmation image is only {:2.0f}% of the "
                  "maximum flux\nin the ACQ/PEAK scan.  Percentages below "
                  "{:2.0f}% often indicate problems\nin the ACQ/PEAK.\n".
                  format(flux_ratio*100, LOW_FLUX_CUTOFF*100))
            badacq |= BAD_RATIO_LOW

        if flux_ratio > HIGH_FLUX_CUTOFF:
            print("The flux in the confirmation image is {:2.0f}% greater than"
                  " the maximum flux\nin the ACQ/PEAK scan.  An excess greater"
                  " than {:3.0f}% indicates\nproblems in the ACQ/PEAK.\n".
                  format((flux_ratio-1)*100, (HIGH_FLUX_CUTOFF-1)*100))
            badacq |= BAD_RATIO_HIGH

        if keywords['goodmax1'] > MAX_GOODMAX:
            badacq |= BAD_SATURATED
            print("Some pixels in the confirmation image were saturated.  "
                  "If saturation also\noccurred in any of the peakup steps, "
                  "it may have affected the centering.\n")

        # Check that the flux level (above pedestal) in the confirmation
        # image is above a minimum value.
        if keywords['obstype'] == "IMAGING":
            if flux < MIN_IMAGING_FLUX:
                print("The flux in the confirmation image is {:2.0f}% of the "
                      "recommended minimum\nof {:.0f} DN for a direct-light "
                      "ACQ/PEAK.  The signal-to-noise in the\nACQ/PEAK may be "
                      "inadequate for an accurate centering.\n".
                      format(flux/MIN_IMAGING_FLUX*100, MIN_IMAGING_FLUX))
                badacq |= BAD_FLUX

        else:
            if flux < MIN_SPECTROSCOPIC_FLUX:
                print("The flux in the confirmation image is {:2.0f}% of the "
                      "recommended minimum\nof {:.0f} DN for a dispersed-light"
                      " ACQ/PEAK.  The signal-to-noise in\nthe ACQ/PEAK may be"
                      " inadequate for an accurate centering.\n".
                      format(flux/MIN_SPECTROSCOPIC_FLUX*100,
                             MIN_SPECTROSCOPIC_FLUX))
                badacq |= BAD_FLUX

        # Search for the word FAILED in keyword PEAKCENT.
        # This will check if flux test failed.
        if "FAILED" in keywords['peakcent']:
            print("The ACQ/PEAK flux test failed, which means that no point in"
                  " the peakup\nscan has a flux that is at least 30% higher "
                  "than any other point.  The\nACQ/PEAK has failed, and the "
                  "telescope has returned to the initial\nposition of the "
                  "ACQ/PEAK\n")
            badacq |= BAD_ACQ

        # If first & last flux values in LINEAR scans are 0.
        if keywords['search'] == "LINEARAXIS1":
            if i_max == 0 or i_max == (keywords['naxis1']-1):
                max_at_end = 1
                badacq |= BAD_END
        elif keywords['search'] == "LINEARAXIS2":
            if j_max == 0 or j_max == (keywords['naxis2']-1):
                max_at_end = 1
                badacq |= BAD_END

        if max_at_end:
            print("The maximum flux in the sequence occurred at one end.\n"
                  "This may indicate that the target was beyond that end\n"
                  "or that a neighboring object affected the acquisition.")

        if badacq == 0:
            print("The confirmation image has a flux between {:3.1f} and "
                  "{:3.1f} times the\nmaximum flux in the peakup, which is "
                  "typical of a successful ACQ/PEAK.".format(LOW_FLUX_CUTOFF,
                                                             HIGH_FLUX_CUTOFF))

    return badacq


def _arcseconds(x):
    """
    Translate pixel values to arcsecond units

    Parameters
    ----------
    x: int, float
        Input pixel value.

    Returns
    -------
    arcsec: float
        Pixel translated to arcsecond unit.
    """

    return x * PLATESCALE


def _v2coord(x, y):
    """
    Translate arcsecond values to v2coordinate system.

    Parameters
    ----------
    x: int, float
        Input arcsecond value.

    y: int, float
        Input arcsecond value.

    Returns
    -------
    v2coord: float
        v2coord value.
    """

    return COSTHETA * x + SINTHETA * y


def _v3coord(x, y):
    """
    Translate arcsecond values to v3coordinate system.

    Parameters
    ----------
    x: int, float
        Input arcsecond value.

    y: int, float
        Input arcsecond value.

    Returns
    -------
    v3coord: float
        v3coord value.
    """

    return COSTHETA * x - SINTHETA * y


def _ndec(x):
    """
    Return input float value rounded up to nearest tenth decimal place.

    Parameters
    ----------
    x: float
        Input value to be rounded up.

    Returns
    -------
    rounded: float
        input value rounded up to nearest tenth decimal place
    """
    if x > 0:
        return int(x*10 + 0.5) / 10.0
    else:
        return int(x*10 - 0.5) / 10.0


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Analyze STIS target acquisition images. :func:`tastis` "
                    "will print general information about each input target "
                    "acquisition image, and will analyze both types of STIS "
                    "target acquisitions: ACQs and ACQ/PEAKs")
    parser.add_argument('filename', type=str,
                        help="Name of the input raw file. For some raw files "
                             "you will need a copy of the spt file in the "
                             "same directory.")
    parser.add_argument('--update', '-u', action='store_true',
                        help='update header')
    args = vars(parser.parse_args())

    tastis(args['filename'], args['update'])
