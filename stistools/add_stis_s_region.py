#!/usr/bin/env python

import glob
import math
import os
import sys
import argparse

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import logging

import pysiaf

__doc__ = """

This script will calculate an S_REGION string for STIS data and assign it
to the S_REGION keyword in the science data header.  If no S_REGION keyword
exists, one will be added to each SCI header, after the existing PA_APER
keyword.  The S_REGION string can be calculated and printed without changing
the science data header by using the --dry_run keyword.

The script uses the pysiaf interface to the HST SIAF database table.  The
PROPAPER keyword and detector are used to look up the entry in the SIAF database
using the correspondence table originally documented in STIS 95-008D ("STIS
Science Apertures Revision D").

The aperture is compared with the extent calculated by using
the science data array and existing WCS and the smaller extent used to calculate
the S_REGION footprint, so that if a subarray was used it will ensure that the
S_REGION footprint doesn't incorrectly specify any sky extent outside the field
of view that is read out.

When used from the command line, the script will calculate and add the S-REGION
keyword to all _raw.fits and _tag.fits files in the current directory.

usage: Add S_REGION value to raw data headers [-h] [--dry_run] [-v]

options:
  -h, --help     show this help message and exit
  --dry_run      Calculate S_REGION value, but don't write to data header[s]

"""

DEGREESTORADIANS = math.pi / 180.0
RADIANTOARCSEC = 180.0 / math.pi * 3600.0

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

STIS_APERTURE_LOOKUP = {'0.1X0.03': '100X030',
                        '0.1X0.06': '100X060',
                        '0.1X0.09': '100X090',
                        '0.1X0.2': '100X200',
                        '0.2X0.05ND': '200X050',
                        '0.2X0.06': '200X060',
                        '0.2X0.06FPA': '200X60A',
                        '0.2X0.06FPB': '200X60B',
                        '0.2X0.06FPC': '200X60C',
                        '0.2X0.06FPD': '200X60D',
                        '0.2X0.06FPE': '200X60E',
                        '0.2X0.09': '200X090',
                        '0.2X0.2': '200X200',
                        '0.2X0.2FPA': '200A',
                        '0.2X0.2FPB': '200B',
                        '0.2X0.2FPC': '200C',
                        '0.2X0.2FPD': '200D',
                        '0.2X0.2FPE': '200E',
                        '0.2X0.5': '200X500',
                        '0.3X0.05ND': '330X050',
                        '0.3X0.06': '330X060',
                        '0.3X0.09': '330X090',
                        '0.3X0.2': '330X200',
                        '0.5X0.5': '500X500',
                        '1X0.06': '1X060',
                        '1X0.2': '1X200',
                        '25MAMA': '25',
                        '25MAMAD1': '25D1',
                        '2X2': '2X2',
                        '31X0.05NDA': '31X050A',
                        '31X0.05NDB': '31X050B',
                        '31X0.05NDC': '31X050C',
                        '36X0.05P45': '36X050P',
                        '36X0.05N45': '36N050N',
                        '36X0.6P45': '36X600P',
                        '36X0.6N45': '36X600N',
                        '50CCD': '50',
                        '50CORON': '50COR',
                        '52X0.05': 'L050',
                        '52X0.05D1': 'L050D1',
                        '52X0.05E1': 'L050E1',
                        '52X0.05F1': 'L050F1',
                        '52X0.05F1-R': 'L050R1',
                        '52X0.05F2': 'L050F2',
                        '52X0.05F2-R': 'L050R2',
                        '52X0.1': 'L100',
                        '52X0.1D1': 'L100D1',
                        '52X0.1E1': 'L100E1',
                        '52X0.1F1': 'L100F1',
                        '52X0.1F1-R': 'L100R1',
                        '52X0.1F2': 'L100F2',
                        '52X0.1F2-R': 'L100R2',
                        '52X0.1B1.0': 'LBAR1',
                        '52X0.1B3.0': 'LBAR3',
                        '52X0.2': 'L200',
                        '52X0.2D1': 'L200D1',
                        '52X0.2E1': 'L200E1',
                        '52X0.2E2': 'L200E2',
                        '52X0.2F1': 'L200F1',
                        '52X0.2F1-R': 'L200R1',
                        '52X0.2F2': 'L200F2',
                        '52X0.2F2-R': 'L200R2',
                        '52X0.5': 'L500',
                        '52X0.5D1': 'L500D1',
                        '52X0.5E1': 'L500E1',
                        '52X0.5E2': 'L500E2',
                        '52X0.5F1': 'L500F1',
                        '52X0.5F1-R': 'L500R1',
                        '52X0.5F2': 'L500F2',
                        '52X0.5F2-R': 'L500R2',
                        '52X2': 'L2',
                        '52X2D1': 'L2D2',
                        '52X2E1': 'L2E1',
                        '52X2E2': 'L2E2',
                        '52X2F1': 'L2F1',
                        '52X2F1-R': 'L2R1',
                        '52X2F2': 'L2F2',
                        '52X2F2-R': 'L2R2',
                        '6X0.06':'6X060',
                        '6X0.2': '6X200',
                        '6X0.5': '6X500',
                        '6X6': '6X6',
                        'BAR5': 'BAR5',
                        'BAR10': 'BAR10',
                        'F25CIII': '25C3',
                        'F25CN182': '25CN182',
                        'F25CN270': '25CN270',
                        'F25LYA': '25LYA',
                        'F25MGII': '25MG2',
                        'F25ND3': '25ND3',
                        'F25ND5': '25ND5',
                        'F25NDQ1': '25NDQ1',
                        'F25NDQ2': '25NDQ2',
                        'F25NDQ3': '25NDQ3',
                        'F25NDQ4': '25NDQ4',
                        'F25QTZ': '25QTZ',
                        'F25QTZD1': '25QTZD1',
                        'F25SRF2': '25SRF2',
                        'F25SRF2D1': '25SRF2D1',
                        'F28X50LP': '28X50LP',
                        'F28X50OII': '28X50O2',
                        'F28X50OIII': '28X50O3',
                        'WEDGEA0.6': 'WGA06',
                        'WEDGEA1.0': 'WGA10',
                        'WEDGEA1.8': 'WGA18',
                        'WEDGEA2.0': 'WGA20',
                        'WEDGEA2.5': 'WGA25',
                        'WEDGEA2.8': 'WGA28',
                        'WEDGEB1.0': 'WGB10',
                        'WEDGEB1.8': 'WGB18',
                        'WEDGEB2.0': 'WGB20',
                        'WEDGEB2.5': 'WGB25',
                        'WEDGEB2.8': 'WGB28'}


def get_files_to_process(rootnames):
    """Create a list of files to process from the list of rootnames

    """

    endings = ['_raw.fits',
               '_tag.fits']

    file_list = []
    for rootname in rootnames:
        if os.path.basename(rootname) != rootname:
            log.warning("{}: rootnames should refer to files in the working directory".format(rootname))
        fitslist = glob.glob(rootname.lower() + '*.fits')
        appended = False
        for input_file in fitslist:
            for ending in endings:
                if input_file.endswith(ending):
                    appended = True
                    file_list.append(input_file)
        if not appended:
            log.warning("No files selected for rootname {}".format(rootname))
    if len(file_list) == 0:
        log.error("No rootnames selected")
    file_list.sort()
    return file_list

def add_s_region(stisfile, hst_siaf, dry_run=False):
    """Calculate the S_REGION keyword for a single STIS file. If the
    dry_run parameter is False, set the S_REGION in the SCI extensions with
    the calculated value.  If keyword isn't present, add it
    """

    open_mode = 'readonly' if dry_run else 'update'
    with fits.open(stisfile, mode=open_mode) as f1:
        log.info('Processing file {}'.format(stisfile))
        hdr0 = f1[0].header
        detector = hdr0['DETECTOR']
        aperture = hdr0['APERTURE']
        propaper = hdr0['PROPAPER']
        siaf_entry = get_siaf_entry(hst_siaf, propaper, detector)
        for ext in f1[1:]:
            if ext.header['EXTNAME'] in ['SCI', 'EVENTS']:
                hdr1 = ext.header
                extname = hdr1['EXTNAME']
                extver = hdr1['EXTVER']
                pa_aper = hdr1['PA_APER']
                pa_aper = pa_aper * DEGREESTORADIANS
                ra_aper = hdr1['RA_APER']
                dec_aper = hdr1['DEC_APER']
                x, y = siaf_entry.closed_polygon_points('idl')
                wcslimits = get_wcs_limits(f1)
                siaflimits = get_siaf_limits(x, y)
                x, y = smallest_size(wcslimits, siaflimits)
                # This is to get the parity right
                x = x * -1.0
                costheta = math.cos(pa_aper)
                sintheta = math.sin(pa_aper)
                dra = x * costheta + y * sintheta
                dra = dra / math.cos(dec_aper * DEGREESTORADIANS) / 3600.0
                ddec = (-x * sintheta + y * costheta) / 3600.0
                ra_corners = ra_aper + dra
                dec_corners = dec_aper + ddec
                s_region = 'POLYGON ICRS'
                for ra, dec in zip(ra_corners, dec_corners):
                    s_region = s_region + ' {} {}'.format(ra, dec)
                log.info("{}[{}, {}] with aperture {} has S_REGION = {}".format(stisfile,
                extname, extver, propaper, s_region))
                if not dry_run:
                    write_keyword_to_header(ext, s_region)
                else:
                    log.info('Dry-run - no changes made to science header')
    return

def write_keyword_to_header(extension, s_region):
    """Write the S_REGION keyword to the header, creating the keyword if it
    doesn't exist"""

    if extension.header.get('S_REGION'):
        # Keyword exists
        extension.header['S_REGION'] = s_region
        log.info('Existing S_REGION keyword updated')
    else:
        # Keyword doesn't exist, need to add it
        extension.header.set('S_REGION', s_region, 'Spatial extent of the observation',
        after='PA_APER')
        log.info('New S_REGION keyword added and updated')
    return

def get_siaf_entry(hst_siaf, aperture, detector):
    """The SIAF entry aperture name doesn't correspond to the STIS
    APERTURE keyword value.  Construct the entry from the aperture and
    detector (here PROPAPER is used as the aperture keyword)

    """
    # All entries start with uppercase o
    entry = 'O'
    # Second letter depends on detector
    detector_letters = {"CCD": "V",
                        "NUV-MAMA": "N",
                        "FUV-MAMA": "F"}
    entry = entry + detector_letters[detector]
    # The rest depends on the aperture
    entry = entry + STIS_APERTURE_LOOKUP[aperture]
    try:
        siaf_entry = hst_siaf[entry]
    except KeyError:
        log.warning("Unable to get SIAF data for entry {}".format(entry))
        log.warning("Trying other wavebands")
        success = False
        for letter in "VNF":
            entry = 'O' + letter + STIS_APERTURE_LOOKUP[aperture]
            try:
                siaf_entry = hst_siaf[entry]
                success = True
                log.info("Succeeded with {}".format(entry))
                break
            except KeyError:
                success = False
        if not success:
            log.error("No matching aperture found")
    return hst_siaf[entry]

def get_wcs_limits(f1):
    # For imaging data, any subarray will sometimes limit the observable
    # footprint.  Calculate the limits of the field of view using the WCS,
    # this will be compared with the limits from projecting the aperture
    hdr1 = f1[1].header
    extname = hdr1['EXTNAME']
    if extname == 'SCI':
        crpix1 = hdr1['CRPIX1']
        crpix2 = hdr1['CRPIX2']
        nx = hdr1['NAXIS1']
        ny = hdr1['NAXIS2']
        if hdr1['CTYPE1'] == 'WAVE':
            xmin = None
            xmax = None
    elif extname == 'EVENTS':
        crpix1 = hdr1['TCRPX2']
        crpix2 = hdr1['TCRPX3']
        nx = hdr1['AXLEN1']
        ny = hdr1['AXLEN2']
        if hdr1['TCTYP2'] == 'WAVE':
            xmin = None
            xmax = None
    else:
        log.error('No SCI or EVENTS extension to get WCS information')
    # Calculate pixel scales in arcsec/pixel in the X (row) and Y (column)
    # direction
    cdelt1, cdelt2 = get_pixel_scales(f1)
    xmin = cdelt1 * crpix1 * -1.0
    xmax = cdelt1 * (nx - crpix1)
    ymin = cdelt2 * crpix2 * -1.0
    ymax = cdelt2 * (ny - crpix2)
    # If the x-axis is dispersed, don't use the WCS limits
    return (xmin, xmax, ymin, ymax)

def get_pixel_scales(f1):
    """Calculate the pixel scales in the X (row) and Y (column) directions
    from the WCS CD matrix, in arcsec/pixel
    """
    hdr = f1[1].header
    extname = hdr['EXTNAME']
    if extname == 'SCI':
        cd1_1 = hdr['CD1_1']
        cd1_2 = hdr['CD1_2']
        cd2_1 = hdr['CD2_1']
        cd2_2 = hdr['CD2_2']
    elif extname == 'EVENTS':
        cd1_1 = hdr['TC2_2']
        cd1_2 = hdr['TC2_3']
        cd2_1 = hdr['TC3_2']
        cd2_2 = hdr['TC3_3']
    else:
        log.error('No SCI or EVENTS extension to get WCS information')
    cdelt1 = math.sqrt(cd1_1*cd1_1 + cd2_1*cd2_1)
    cdelt2 = math.sqrt(cd1_2*cd1_2 + cd2_2*cd2_2)
    cdelt1 = cdelt1 * 3600.0
    cdelt2 = cdelt2 * 3600.0
    return cdelt1, cdelt2

def get_siaf_limits(x, y):
    """Given arrays of X and Y values from the SIAF in ideal coordinates
    (relative to the aperture reference point), return the X and Y limits
    """
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    return (xmin, xmax, ymin, ymax)

def smallest_size(wcslimits, siaflimits):
    """Find the limiting values of X and Y (ideal coords) given
    the limits from the SIAF projected aperture and those calculated from the
    WCS applied to the limits of the data array
    """
    if wcslimits[0] == None:
        xmin = siaflimits[0]
        xmax = siaflimits[1]
    else:
        xmin = max(wcslimits[0], siaflimits[0])
        xmax = min(wcslimits[1], siaflimits[1])
    ymin = max(wcslimits[2], siaflimits[2])
    ymax = min(wcslimits[3], siaflimits[3])
    return (np.array([xmin, xmin, xmax, xmax, xmin]),
            np.array([ymin, ymax, ymax, ymin, ymin]))

def coords_from_s_region(s_region):
    """Helper function to extract the RA and DEC as lists of floating
    point numbers from the S_REGION string"""

    first_coordinate = 0
    coords = s_region.split()
    while coords[first_coordinate] in ['POLYGON', 'ICRS']:
        first_coordinate = first_coordinate + 1
    ra = coords[first_coordinate::2]
    ra = [float(x) for x in ra]
    dec = coords[first_coordinate+1::2]
    dec = [float(x) for x in dec]
    return (ra, dec)

def write_region_file(input_file, include_mast=True):
    """Convenience function to write a DS9 regions file from the S_REGION
    value in the header of a file.  If include_mast is set to True, will
    also get the existing S_REGION value from the mast catalog and include
    it in the region file.  The region file will include the keyword S_REGION
    in blue, the RA_APER and DEC_APER in red, and the mast catalog value of the
    S_REGION in green, if include_mast is True"""

    f1 = fits.open(input_file)
    rootname = f1[0].header['ROOTNAME']
    ra_aper = f1[1].header['RA_APER']
    dec_aper = f1[1].header['DEC_APER']
    s_region = f1[1].header['S_REGION']
    f1.close()
    f_region = open(rootname+'.reg', mode='w')
    text = 'icrs; polygon'
    ra, dec = coords_from_s_region(s_region)
    for x, y in zip(ra, dec):
        text = text + ' {} {}'.format(x, y)
    text = text + " # color=blue"
    f_region.write(text)
    f_region.write('\n')
    text = 'icrs; point' + ' {} {}'.format(ra_aper, dec_aper) + ' # color=red \n'
    # Get existing S_REGION value from MAST
    if include_mast:
        from astroquery.mast import Observations
        row = Observations.query_criteria(obs_id=rootname)
        existing_s_region = row['s_region'][0]
        f_region.write(text)
        text = 'icrs; polygon'
        ra, dec = coords_from_s_region(existing_s_region)
        for x, y in zip(ra, dec):
            text = text + ' {} {}'.format(x, y)
        text = text + " # color=green"
        f_region.write(text)
        f_region.write('\n')
    f_region.close()
    return



def main(rootnames, dry_run=False):
    if rootnames is None:
        log.error("No rootnames specified")
        return
    files_to_process = get_files_to_process(rootnames)
    hst_siaf = pysiaf.Siaf('HST')
    for input_file in files_to_process:
        add_s_region(input_file, hst_siaf, dry_run=dry_run)
    return

def call_main():

    parser = argparse.ArgumentParser(
    """Add S_REGION value to raw data headers"""
    )

    parser.add_argument('rootnames', nargs='+',
        help='Rootnames to be processed')

    parser.add_argument(
        '--dry_run', action='store_true',
        help="Calculate S_REGION value, but don't write to data header[s]")

    args = parser.parse_args()
    # if '--version' in args:
    #     print(__version__)
    #     sys.exit(0)

    main(args.rootnames, dry_run=args.dry_run)

if __name__ == '__main__':
    call_main()
