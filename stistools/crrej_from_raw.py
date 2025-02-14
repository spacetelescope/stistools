#!/usr/bin/env python

import numpy as np
from tempfile import TemporaryDirectory
import os
import sys
import warnings
from numbers import Number
from astropy.io import fits
from astropy.table import Table
from .calstis import calstis
from .r_util import expandFileName


__taskname__ = "crrej_from_raw"
__version__ = "0.1"
__vdate__ = "12-May-2023"
__author__ = "Sean Lockwood, STScI, May 2023."

VERBOSE = False
HIST_LINE = 'CRREJTAB temporarily modified for processing.'
ENV_VAR = 'TEMPCRR'


def crrej_from_raw(input, wavecal='', outroot='', savetmp=False, verbose=False,
                   timestamps=False, trailer='', print_version=False, print_revision=False,
                   crrejtab='', scalense=None, initgues=None, skysub=None, crsigmas=None,
                   crradius=None, crthresh=None, badinpdq=None, crmask=None):
    '''Calibrate STIS data (calstis) with modified cosmic ray rejection parameters.

    A temporary CRREJTAB is created, with values derived from user-supplied inputs and
    defaults derived from the file's ext=0 CRREJTAB keyword.  This default file may
    be overridden by specifying a different `crrejtab`.

    This routine should be used instead of `stistools.ocrreject.ocrreject()` when the
    order of operation is of concern.

    **Note:**
    This routine does not replace `ocrreject`'s functionality for combining multiple
    CCD datasets together into one CRJ file (i.e. the "all" parameter).

    Parameters
    ----------
    input: str
        Name of the input file.

    wavecal: str, default=''
        Name of the input wavecal file, or '' (the default).  This is
        only needed if the name is not the "normal" name
        (rootname_wav.fits).

    outroot: str, default=''
        Root name for the output files, or '' (the default).  This can
        be a directory name, in which case the string must end in '/'.

    savetmp: bool, default=False
        True if calstis should not delete temporary files.

    verbose: bool, default=False
        If True, calstis will print more info.

    timestamps: bool, default=False
        If True, calstis will print the date and time at various points
        during processing.

    trailer: str, default=''
        If specified, the standard output and standard error will be
        written to this file instead of to the terminal.  Note, however,
        that if print_version or print_revision is specified, the value
        will be printed to the terminal, and any name given for the
        trailer will be ignored.

    print_version: bool, default=False
        If True, calstis will print the version number (a string) and
        then return 0.

    print_revision: bool, default=False
        If True, calstis will print the full version string and then
        return 0.

    crrejtab: str, default=''
        This argument may be used to override the CRREJTAB value in the
        primary header of the input files.  This file determines the
        default CR-rejection parameters below.

    scalense: str, default=''
        Multiplicative scale factor applied to noise.
        If specified, this overrides SCALENSE in the CRREJTAB.

    initgues: str, default=''
        Initial guess method.
        If specified, this overrides INITGUES in the CRREJTAB.
        The allowed values are 'minimum' and 'median' and ''.

    skysub: str, default=''
        Sky value subtracted.
        If specified, this overrides SKYSUB in the CRREJTAB.
        The allowed values are 'none', 'mode' and ''.

    crsigmas: str, default=''
        Statistical rejection criteria.
        If specified, this overrides CRSIGMAS in the CRREJTAB.  The
        value should be a comma-separated string of one or more
        integer or float values.  For each such value, calstis will
        perform one cosmic-ray-rejection cycle, with the sigma taken
        from the numerical value that was specified.

    crradius: float or None, default=None
        If not None, this overrides CRRADIUS in the CRREJTAB.  This is
        the rejection propagation radius in pixels (e.g. 1.5).  After
        finding an outlier (a cosmic ray hit), adjacent pixels can also
        be flagged and excluded.  Neighboring pixels will be rejected if
        their values are discrepant by more than crthresh * sigmas * noise,
        where noise is based on the noise model (i.e. Poisson noise and
        readout noise).

    crthresh: float or None, default=None
        If not None, this overrides CRTHRESH in the CRREJTAB.  This is the
        rejection propagation threshold (e.g. 0.8).  If crthresh = 0 then
        all adjacent pixels (see crradius) will be rejected.

    badinpdq: int or None, default=None
        If specified, this overrides BADINPDQ in the CRREJTAB.  This is a
        data quality flag (or bitwise OR of flags) to allow rejection of
        pixels in the input images when forming the "guess" image (the
        image with which to compare the input images when looking for
        outliers).

    crmask: bool or None, default=None
        If specified, this overrides CRMASK in the CRREJTAB.  crmask = True
        means that the cosmic rays that are detected should be flagged in
        the DQ (data quality) extensions of the input files.

    Returns
    -------
    status: int
        0 is OK.
        1 is returned if cs0.e (the calstis host executable) returned a
        non-zero status.  If verbose is True, the value returned by cs0.e
        will be printed.
        2 is returned if the specified input file or files were not found.
    '''
    global VERBOSE
    VERBOSE = verbose

    if print_version or print_revision:
        calstis(None, print_version=print_version, print_revision=print_revision)
        sys.exit(0)

    crr = determine_crrejtab(input, crrejtab=crrejtab)

    # Determine CRSPLIT & MEANEXP to determine appropriate row in CRR:
    with fits.open(input) as f:
        exts = range(1, len(f), 3)
        crsplit = len(exts)
        totexptime = 0.
        for ext in exts:
            totexptime += f[ext].header.get('EXPTIME', np.nan)

    meanexp = totexptime / crsplit
    if np.isnan(meanexp):
        raise ValueError('Problem finding EXPTIME keywords in SCI extensions of input file.')

    # Locate row corresponding to the input RAW file:
    crr = crr[(crr['CRSPLIT'] == crsplit) & (crr['MEANEXP'] >= meanexp)]
    crr = crr[np.argmin(crr['MEANEXP'] - meanexp)]

    crr_par = {}  # Final parameters that will be used to generate a temporary CRREJ table
    crr_par['crsplit'] = crr['CRSPLIT']  # Used for row selection in calstis
    crr_par['meanexp'] = crr['MEANEXP']  # Used for row selection in calstis

    # Determine either user-specified values or defaults from CRR file:
    # Checking for Number type allows user to specify 0 or 0.0.
    crr_par['scalense'] = str(scalense if isinstance(scalense, Number)
                              else (scalense or crr['SCALENSE'])).strip()
    crr_par['initgues'] = str(initgues or crr['INITGUES']).lower().strip()
    crr_par['skysub'] = str(skysub or crr['SKYSUB']).lower().strip()
    crr_par['crsigmas'] = str(float(crsigmas) if isinstance(crsigmas, Number)
                              else (crsigmas or crr['CRSIGMAS'])).strip()
    crr_par['crradius'] = float(crradius if isinstance(crradius, Number)
                                else (crradius or crr['CRRADIUS']))
    crr_par['crthresh'] = float(crthresh if isinstance(crthresh, Number)
                                else (crthresh or crr['CRTHRESH']))
    crr_par['badinpdq'] = int(badinpdq if isinstance(badinpdq, Number)
                              else (badinpdq or crr['BADINPDQ']))
    crr_par['crmask'] = bool(crmask if (crmask is not None) else crr['CRMASK'])

    # Check some input constraints:
    if crr_par['initgues'] not in {'minimum', 'median'}:
        raise ValueError(f'INITGUES="{crr_par["initgues"]}" not a permitted '
                         'value of "minimum" or "median".')
    if crr_par['skysub'] not in {'mode', 'none'}:
        raise ValueError(f'SKYSUB="{crr_par["skysub"]}" not a permitted '
                         'value of "mode" or "none".')

    # Convert dict to Astropy Table:
    new_crr = create_new_crr(crr_par)

    with TemporaryDirectory(prefix='crrej_from_raw_') as directory:
        # Write new CRREJTAB to a temporary location:
        new_crr_name = 'temp_crr.fits'
        new_crr.write(os.path.join(directory, new_crr_name))
        if VERBOSE:
            print(new_crr)

        # Modify the input file to use this temporary CRREJTAB file:
        os.environ[ENV_VAR] = directory  # + os.path.sep
        with fits.open(input, 'update') as f:
            if not f[0].header.get('CRREJTAB', 'N/A').startswith('$' + ENV_VAR):
                f[0].header.set('PREV_CRR', after='CRREJTAB',
                                value=f[0].header['CRREJTAB'],
                                comment='previous CRREJTAB')
            f[0].header['CRREJTAB'] = f'${ENV_VAR}/{new_crr_name}'

        try:
            # Call calstis with remaining input parameters:
            status = calstis(input, wavecal=wavecal, outroot=outroot, savetmp=savetmp,
                             verbose=verbose, timestamps=timestamps, trailer=trailer,
                             print_version=print_version, print_revision=print_revision)
        finally:
            # Clean up input RAW file:
            with fits.open(input, 'update') as f:
                f[0].header['CRREJTAB'] = f[0].header['PREV_CRR']
                del f[0].header['PREV_CRR']

                # Revert HISTORY in RAW file:
                try:
                    # Find last instance of HIST_LINE in the HISTORY:
                    r = [i for i, x in enumerate(f[0].header['HISTORY'])
                         if x.strip() == HIST_LINE][-1]
                    del f[0].header['HISTORY', r]
                except (IndexError, KeyError):
                    pass

    return status


def determine_crrejtab(input, crrejtab=None):
    '''Parse the input file's ext=0 header for the CRREJTAB.
    A user-specified CRREJTAB may be provided instead.

    - update the input file's HISTORY
    - update the input file's CRREJTAB keyword to point to a new temporary file
    - temporarily store the old CRREJTAB value (if valid) in a new PREV_CRR keyword
    - read the CRREJTAB FITS table

    Parameters
    ----------
    input: str
        Name of the RAW FITS file to update.

    crrejtab: str or None, default=None
        User-supplied CRREJTAB to be used instead of the file specified in the ext=0
        header.

    Returns
    -------
        Astropy Table of CRREJTAB contents
    '''
    with fits.open(input, 'update') as f:
        # Add a HISTORY line to the RAW file so subsequent products get modified.
        # This line will be removed from the RAW file later when it is reverted.
        f[0].header.add_history(HIST_LINE)

        if not crrejtab:
            # Determine CRREJTAB in order to find applicable row and determine
            # default parameters:
            crrejtab = f[0].header.get('CRREJTAB', 'N/A')
            if crrejtab.startswith('$' + ENV_VAR) and ('PREV_CRR' in f[0].header):
                # Typically CRREJTAB will revert and PREV_CRR will get deleted, but
                # if the code exits early the RAW file may need to be repaired.
                if VERBOSE:
                    print(f'Using PREV_CRR="{f[0].header["PREV_CRR"]}" for CRREJTAB...')
                crrejtab = f[0].header['CRREJTAB'] = f[0].header['PREV_CRR']

    crrejtab = expandFileName(crrejtab)  # Works with both "$dir/fn" and "oref$fn" syntax
    if crrejtab.strip().upper() == 'N/A':
        raise ValueError('CRREJTAB not specified.')
    if not os.access(crrejtab, os.F_OK):
        raise FileNotFoundError(f'CRREJTAB ({crrejtab}) not found.')
    with warnings.catch_warnings(record=True) as _:
        # Silence warnings related to malformed FITS keywords in the CRDS CRREJTAB.
        crr = Table.read(crrejtab, hdu=1)

    return crr


def create_new_crr(crr_par):
    '''Convert a dict to a 1-row Astropy Table with the correct column dtypes.

    Parameters
    ----------
    crr_par: dict
        Dictionary of CRREJECTAB contents

    Returns
    -------
    Astropy Table of CRREJECTAB contents
    '''
    # Convert to an Astropy table:
    type_map = {'CRSPLIT': 'i2', 'MEANEXP': 'f', 'SCALENSE': 'a8', 'INITGUES': 'a8',
                'SKYSUB': 'a4', 'CRSIGMAS': 'a20', 'CRRADIUS': 'f', 'CRTHRESH': 'f',
                'BADINPDQ': 'l', 'CRMASK': bool, }

    if {x.upper() for x in crr_par} != {y.upper() for y in type_map}:
        raise ValueError('New CRREJTAB not specified correctly.')

    dtypes = [type_map[k.upper()] for k in crr_par]
    new_crr = Table({k.upper(): [v] for k, v in crr_par.items()}, dtype=dtypes)

    return new_crr
