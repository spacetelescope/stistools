"""HSTCAL regression test helpers."""
from six.moves import urllib
import getpass
import os
import sys
import math
from io import StringIO
import shutil
import datetime
from os.path import splitext
from difflib import unified_diff

import pytest
import requests
from astropy.io import fits
from astropy.io.fits import FITSDiff
from astropy.table import Table
from astropy.utils.data import conf


from .helpers.io import get_bigdata, upload_results

__all__ = ['download_crds',
           'ref_from_image', 'raw_from_asn', 'BaseACS',
           'BaseSTIS', 'BaseWFC3IR', 'BaseWFC3UVIS', 'BaseWFPC2']


def _download_file(url, filename, filemode='wb', timeout=None):
    """Generic remote data download."""
    if url.startswith('http'):
        r = requests.get(url, timeout=timeout)
        with open(filename, filemode) as fout:
            fout.write(r.content)
    elif url.startswith('ftp'):  # TODO: Support filemode and timeout.
        urllib.request.urlretrieve(url, filename=filename)
    else:  # pragma: no cover
        raise ValueError('Unsupported protocol for {}'.format(url))


def download_crds(refdir, refname, timeout=None):
    """Download a CRDS file from HTTP or FTP to current directory."""
    # CRDS file for given name never changes, so no need to re-download.
    if os.path.exists(refname):
        return

    try:
        url = 'http://ssb.stsci.edu/cdbs/{}/{}'.format(refdir, refname)
        local_file = os.path.abspath(refname)
        print("Downloading CRDS file: {}".format(local_file))
        _download_file(url, refname, timeout=timeout)
    except Exception:  # Fall back to FTP
        url = 'ftp://ftp.stsci.edu/cdbs/{}/{}'.format(refdir, refname)
        _download_file(url, refname, timeout=timeout)


def _get_reffile(hdr, key):
    """Get ref file from given key in given FITS header."""
    ref_file = None
    if key in hdr:  # Keyword might not exist
        ref_file = hdr[key].strip()
        if ref_file.upper() == 'N/A':  # Not all ref file is defined
            ref_file = None
    return ref_file


def ref_from_image(input_image):
    """
    Return a list of reference filenames, as defined in the primary
    header of the given input image, necessary for calibration; i.e.,
    only those associated with ``*CORR`` set to ``PERFORM`` will be
    considered.
    """
    # NOTE: Add additional mapping as needed.
    # Map mandatory CRDS reference file for instrument/detector combo.

    reffile_lookup = ['BPIXTAB', 'DARKFILE', 'PFLTFILE', 'LFLTFILE', 'PHOTTAB',
                      'IMPHTTAB', 'APERTAB', 'CCDTAB', 'BIASFILE', 'CRREJTAB',
                      'IDCTAB', 'TDSTAB', 'SPTRCTAB', 'SDCTAB', 'PHOTTAB',
                      'PCTAB', 'TDCTAB', 'MLINTAB', 'GACTAB', 'WCPTAB',
                      'LAMPTAB', 'APDESTAB', 'XTRACTAB', 'DISPTAB', 'INANGTAB',
                      'CDSTAB', 'ECHSCTAB', 'EXSTAB', 'HALOTAB', 'TELTAB',
                      'RIPTAB', 'SRWTAB']

    ref_files = []
    hdr = fits.getheader(input_image, ext=0)

    for reffile in reffile_lookup:
        s = _get_reffile(hdr, reffile)
        if s is not None:
            ref_files.append(s)

    return ref_files


def raw_from_asn(asn_file, suffix='_raw.fits'):
    """Return a list of RAW input files in a given ASN."""
    raw_files = []
    tab = Table.read(asn_file, format='fits')

    for row in tab:
        if row['MEMTYPE'].startswith('PROD'):
            continue
        pfx = row['MEMNAME'].lower().strip().replace('\x00', '')
        raw_files.append(pfx + suffix)

    return raw_files


# Base classes for actual tests.
# NOTE: Named in a way so pytest will not pick them up here.
#@pytest.mark.bigdata
class BaseCal(object):
    prevdir = os.getcwd()
    use_ftp_crds = True
    timeout = 30  # seconds
    tree = ''
    results_root = 'datb-stistools/results'

    # Numpy default for allclose comparison
    rtol = 5e-7

    atol = 0

    # To be defined by instrument
    refstr = ''
    prevref = ''
    input_loc = ''
    ref_loc = ''
    ignore_keywords = []

    # To be defined by individual test
    subdir = ''

    @pytest.fixture(autouse=True)
    def setup_class(self, tmpdir, envopt):
        """
        Run test in own dir so we can keep results separate from
        other tests.
        """
        if not tmpdir.ensure(self.subdir, dir=True):
            p = tmpdir.mkdir(self.subdir).strpath
        else:
            p = tmpdir.join(self.subdir).strpath
        os.chdir(p)

        # NOTE: This could be explicitly controlled using pytest fixture
        #       but too many ways to do the same thing would be confusing.
        #       Refine this logic if using pytest fixture.
        # HSTCAL cannot open remote CRDS on FTP but central storage is okay.
        # So use central storage if available to avoid FTP.
        if self.prevref is None or self.prevref.startswith(('ftp', 'http')):
            os.environ[self.refstr] = p + os.sep
            self.use_ftp_crds = True

        # Turn off Astrometry updates
        os.environ['ASTROMETRY_STEP_CONTROL'] = 'OFF'

        # This controls astropy.io.fits timeout
        conf.remote_timeout = self.timeout

        # Update tree to point to correct environment
        self.tree = envopt

    def teardown_class(self):
        """Reset path and variables."""
        conf.reset('remote_timeout')
        os.chdir(self.prevdir)
        if self.use_ftp_crds and self.prevref is not None:
            os.environ[self.refstr] = self.prevref

    def get_data(self, *args):
        """
        Download `filename` into working directory using
        `helpers/io/get_bigdata`.  This will then return the full path to
        the local copy of the file.
        """
        local_file = get_bigdata(self.tree, self.input_loc, *args)

        return local_file

    def get_input_file(self, *args, refsep='$'):
        """
        Download or copy input file (e.g., RAW) into the working directory.
        The associated CRDS reference files in ``refstr`` are also
        downloaded, if necessary.
        """
        filename = self.get_data(*args)

        print(filename)

        ref_files = ref_from_image(filename)
        print("Looking for REF_FILES: {}".format(ref_files))

        for ref_file in ref_files:
            if ref_file.strip() == '':
                continue
            if refsep not in ref_file:  # Local file
                refname = self.get_data('customRef', ref_file)
            else:  # Download from FTP, if applicable
                s = ref_file.split(refsep)
                refdir = s[0]
                refname = s[1]
                if self.use_ftp_crds:
                    download_crds(refdir, refname, timeout=self.timeout)
        return filename

    def compare_outputs(self, outputs, raise_error=True, delete_history=False):
        """
        Compare output with "truth" using appropriate
        diff routine; namely,
            ``fitsdiff`` for FITS file comparisons
            ``unified_diff`` for ASCII products.
        Parameters
        ----------
        outputs : list of tuple
            A list of tuples, each containing filename (without path)
            of CALXXX output and truth, in that order.
        raise_error : bool
            Raise ``AssertionError`` if difference is found.
        Returns
        -------
        report : str
            Report from ``fitsdiff``.
            This is part of error message if ``raise_error=True``.
        """
        all_okay = True
        creature_report = ''
        # Create instructions for uploading results to artifactory for use
        # as new comparison/truth files
        testpath, testname = os.path.split(os.path.abspath(os.curdir))
        # organize results by day test was run...could replace with git-hash
        whoami = getpass.getuser() or 'nobody'
        dt = datetime.datetime.now().strftime("%d%b%YT")
        ttime = datetime.datetime.now().strftime("%H_%M_%S")
        user_tag = 'NOT_CI_{}_{}'.format(whoami, ttime)
        build_tag = os.environ.get('BUILD_TAG',  user_tag)
        build_suffix = os.environ.get('BUILD_MATRIX_SUFFIX', 'standalone')
        testdir = "{}_{}_{}".format(testname, build_tag, build_suffix)
        tree = os.path.join(self.results_root, self.input_loc,
                            dt, testdir) + os.sep

        updated_outputs = []
        for actual, desired in outputs:
            # Get "truth" image
            s = self.get_data('truth', desired)
            if s is not None:
                desired = s

            if actual.endswith('fits'):
                # Working with FITS files...
                if delete_history is True:
                    actual = fits.open(actual)
                    desired = fits.open(desired)
                    if 'HISTORY' in actual[0].header:
                        del actual[0].header['HISTORY']
                    if 'HISTORY' in desired[0].header:
                        del desired[0].header['HISTORY']

                fdiff = FITSDiff(actual, desired, rtol=self.rtol,
                                 atol=self.atol,
                                 ignore_keywords=self.ignore_keywords)

                if delete_history is True:
                    actual.close()
                    desired.close()

                creature_report += fdiff.report()
                if not fdiff.identical:
                    # Only keep track of failed results which need to
                    # be used to replace the truth files (if OK).
                    updated_outputs.append((actual, desired))
                if not fdiff.identical and all_okay:
                    all_okay = False
            else:
                # Try ASCII-based diff
                with open(actual) as afile:
                    actual_lines = afile.readlines()
                with open(desired) as dfile:
                    desired_lines = dfile.readlines()
                udiff = unified_diff(actual_lines, desired_lines,
                                     fromfile=actual, tofile=desired)

                old_stdout = sys.stdout
                udiffIO = StringIO()
                sys.stdout = udiffIO
                sys.stdout.writelines(udiff)
                sys.stdout = old_stdout
                udiff_report = udiffIO.getvalue()
                creature_report += udiff_report
                if len(udiff_report) > 2 and all_okay:
                    all_okay = False
                if len(udiff_report) > 2:
                    # Only keep track of failed results which need to
                    # be used to replace the truth files (if OK).
                    updated_outputs.append((actual, desired))

        if not all_okay:
            # Write out JSON file to enable retention of different results
            new_truths = [os.path.abspath(i[1]) for i in updated_outputs]
            for files in updated_outputs:
                print("Renaming {} as new 'truth' file: {}".format(
                      files[0], files[1]))
                shutil.move(files[0], files[1])
            log_pattern = [os.path.join(os.path.dirname(x), '*.log')
                           for x in new_truths]
            upload_results(pattern=new_truths + log_pattern,
                           testname=testname,
                           target=tree)

        if not all_okay and raise_error:
            raise AssertionError(os.linesep + creature_report)

        return creature_report


class BaseSTIS(BaseCal):

    refstr = 'oref'
    prevref = os.environ.get(refstr)
    input_loc = ''
    ref_loc = '/ref'
    ignore_keywords = ['date', 'filename', 'iraf-tlm', 'fitsdate', 'history']
                #''cal_ver']

    def read_image(self, filename):
        """
        Read the image from a fits file
        """
        hdu = fits.open(filename)

        image = hdu[1].data
        hdu.close()
        return image


def add_suffix(fname, suffix, range=None):
    """Add suffix to file name
    Parameters
    ----------
    fname: str
        The file name to add the suffix to
    suffix: str
        The suffix to add_suffix
    range: range
        If specified, the set of indexes will be added to the
        outputs.
    Returns
    -------
    fname, fname_with_suffix
        2-tuple of the original file name and name with suffix.
        If `range` is defined, `fname_with_suffix` will be a list.
    """
    fname_root, fname_ext = splitext(fname)
    if range is None:
        with_suffix = ''.join([
            fname_root,
            '_',
            suffix,
            fname_ext
        ])
    else:
        with_suffix = []
        for idx in range:
            with_suffix.append(''.join([
                fname_root,
                '_',
                str(idx),
                '_',
                suffix,
                fname_ext
            ]))

    return fname, with_suffix
