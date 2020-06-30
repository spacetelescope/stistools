#! /usr/bin/env python

import os
import shutil
import stat
import re
from astropy.io import fits
from tempfile import mkdtemp
import warnings

from ..r_util import expandFileName
from ..calstis import calstis


def prepspec(inspec, outroot='./', darkfile=None, pixelflat=None, initguess=None):
    """Calibrate STIS CCD G750L or G750M spectrum before defringing.

    Based on the PyRAF `stsdas.hst_calib.stis.prepspec` task.

    Parameters
    ----------
    inspec: str
        Name of input 'raw' science spectrum
    outroot: str
        Root for output file name.  (Default='./')
    darkfile: str or None
        Name of superdark image.  If None, use DARKFILE in main header of input spectrum.
    pixelflat: str or None
        Name of pixel-to-pixel flat.  If None, use PIXELFLAT in main header of input
        spectrum.
    initguess: str or None
        Method for initial value estimate for `ocrreject`: {None, 'minimum', 'median'}.
        (Default=None; Use the value in the CRREJTAB.)

    Returns
    -------
    outname: str
        Fully qualified name of prepared spectrum (CRJ or SX2 file)

    """
    # These notes are based on the old STSDAS algorithm.  It may make sense to change the 
    # order or processing when converting to Python.

    # Note that stistools.basic2d.basic2d(), stistools.ocrreject.ocrreject(), and 
    # stistools.x2d.x2d() likely ignore ext=0 calibration flags.

    # Check inputs:
    #    - inspec:  HST/STIS, G750M/G750L, correct filetype, not already corrected
    #    - outroot:  handle dir vs filename; ".fits" in name or not
    #    - don't allow inspec to be of type {crj, sx2}
    #    - print/log name of resolved reference files
    #    - read in reference file data

    # Check if we can we CR-combine the data:
    #    - Check inspec has more than 1 {SCI, ERR, DQ} ext groups (NEXTEND//3).
    #    - Handle CR-rejection header keywords when NRPTEXP != NEXTEND//3.

    # perform only {DQICORR, BLEVCORR} via stistools.basic2d.basic2d() --> Produces temporary file
    # perform CRCORR via stistools.ocrreject.ocrreject() --> Produces CRJ file
    # perform only {BIASCORR, DARKCORR, FLATCORR, PHOTCORR; STATFLAG} via stistools.basic2d.basic2d()
    # If G750L:
    #    return CRJ file
    # If G750M:
    #    perform only {HELCORR, X2DCORR} via stistools.x2d.x2d()
    #    return SX2 file

    science_data = os.path.abspath(expandFileName(inspec))
    sci_root = re.split('\.fits.*', os.path.basename(science_data),
        flags=re.IGNORECASE)[0].rsplit('_', 1)[0]
    opt_elem = fits.getval(science_data, 'OPT_ELEM').strip().upper()
    if (fits.getval(science_data, ext=0, keyword='DETECTOR').strip().upper() != 'CCD') or \
            (fits.getval(science_data, ext=0, keyword='INSTRUME').strip() != 'STIS') or \
            (opt_elem not in ['G750L', 'G750M']):
        raise ValueError('prepspec:  Intended for use on STIS/CCD G750L & G750M data!')
    if (initguess is not None) and (initguess.lower() not in ['minimum', 'median']):
        raise ValueError('initguess must be in {None, "minimum", "median"}!')
    if os.path.isdir(outroot):
        outroot = os.path.normpath(outroot) + os.sep
    if not os.access(os.path.dirname(outroot), os.W_OK):
        raise IOError('Cannot write to:  {}'.format(os.path.dirname(outroot)))
    #if (not os.path.isdir(outroot)) and os.access(outroot, os.F_OK):
    #    raise FileExistsError('Previous outroot already exists:  {}'.format(outroot))

    # Make sure the necessary header keywords are set to PERFORM:
    with fits.open(science_data, 'update') as f:
        f[0].header['STATFLAG'] = True
        for keyword in ['DQICORR', 'BLEVCORR', 'BIASCORR', 'DARKCORR', 'FLATCORR', 'CRCORR']:
            if not f[0].header[keyword].upper().startswith('COMPLETE'):
                f[0].header[keyword] = 'PERFORM'

        # A few extra calibration steps if G750M data:
        if opt_elem == 'G750M':
            for keyword in ['HELCORR', 'X2DCORR']:
                if not f[0].header[keyword].upper().startswith('COMPLETE'):
                    f[0].header[keyword] = 'PERFORM'
            # For comparison with PyRAF/IRAF products:
            if f[0].header['WAVECORR'].upper().startswith('PERFORM'):
                f[0].header['WAVECORR'] = 'OMIT'

        ref_types = {
            'DARKFILE': os.path.abspath(darkfile  or expandFileName(f[0].header['DARKFILE'])),
            'PFLTFILE': os.path.abspath(pixelflat or expandFileName(f[0].header['PFLTFILE'])),}

        # Handle non-default CR-rejection initial guesses:
        if initguess:
            # Copy the old CRREJTAB to a temporary location:
            crrejtab_tmp_dir = mkdtemp(prefix='crrejtab_')
            orig_crrejtab_name = expandFileName(f[0].header['CRREJTAB'])
            try:
                shutil.copy(orig_crrejtab_name, crrejtab_tmp_dir)
            except FileNotFoundError:
                warnings.warn('\nUsing CRREJTAB in $oref instead of ext=0 CRREJTAB="{}"'.format(f[0].header['CRREJTAB']))
                orig_crrejtab_name = os.path.join(os.environ['oref'], os.path.basename(orig_crrejtab_name))
                shutil.copy(orig_crrejtab_name, crrejtab_tmp_dir)
            new_crrejtab_name = os.path.join(crrejtab_tmp_dir, os.path.basename(orig_crrejtab_name))
            os.chmod(new_crrejtab_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

            # Update header of science_data file:
            ref_types['CRREJTAB'] = os.path.abspath(new_crrejtab_name)

            # Modify the temporary CRREJTAB:
            with fits.open(new_crrejtab_name, 'update') as crrejtab:
                crrejtab[1].data['INITGUES'] = initguess.lower()
        else:
            crrejtab_tmp_dir = None

        # Populate/repopulate the inflat header accordingly:
        for i, (ref_type, ref) in enumerate(ref_types.items()):
            if not os.access(ref, os.F_OK):
                raise FileNotFoundError('Cannot access reference file:  {}'.format(ref))

            # Handle reference file paths via environment variables:
            ref_var = 'reff{:.0f}'.format(i+1)
            os.environ[ref_var] = os.path.abspath(os.path.dirname(ref)) + os.sep
            # Keep $oref where it's the same:
            if os.path.normpath(os.environ[ref_var]) == os.path.normpath(os.environ['oref']):
                ref_var = 'oref'
                f[0].header[ref_type] = '{}${}'.format(ref_var, os.path.basename(ref))
            else:
                f[0].header[ref_type] = '${}/{}'.format(ref_var, os.path.basename(ref))

    # Calibrate with calstis:
    trl_file = '{}_trl.txt'.format(sci_root)  # Output log goes here
    cwd = os.getcwd()
    try:
        # Run calstis from within the directory with the data to find EPC files properly:
        os.chdir(os.path.dirname(science_data))
        res = calstis(os.path.basename(science_data), outroot=outroot, trailer=trl_file)
    finally:
        os.chdir(cwd)
        if crrejtab_tmp_dir:
            shutil.rmtree(crrejtab_tmp_dir)  # Temporary CRREJTAB

    # Print out calstis log:
    trl_file = os.path.join(os.path.dirname(science_data), trl_file)
    with open(trl_file) as trl:
        for line in trl:
            print ('    ' + line.rstrip())

    # Raise exception on bad calstis exit code:
    if res != 0:
        raise Exception('CalSTIS exited with code {}'.format(res))


def call_prepspec():
    """Command line entry point for prepspec().
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Correct STIS CCD G750L or G750M spectrum for fringing')
    parser.add_argument('inspec', type=str, help='Name of input "raw" science spectrum')
    parser.add_argument('outroot', type=str, default='./',
        help='Root for output file name. (Default="./")')
    parser.add_argument('--darkfile', '-d', type=str, default=None,
        help='Name of superdark image. If omitted, use DARKFILE in main header of input '
             'spectrum.')
    parser.add_argument('--pixelflat', '-f', type=str, default=None,
        help='Name of pixel-to-pixel flat.  If omitted, use PIXELFLAT in main header of '
             'input spectrum.')
    parser.add_argument('--initguess', '-i', type=str, default=None,
        choices=['None', 'minimum', 'median'],
        help='Method for initial value estimate for ocrreject (Default=None; use value from CRREJTAB)')

    args = vars(parser.parse_args())
    if isinstance(args['initguess'], str) and args['initguess'].lower() == 'none':
        args['initguess'] = None
    prepspec(**args)


if __name__ == '__main__':
    call_prepspec()
