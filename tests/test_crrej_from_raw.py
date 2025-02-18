import pytest
import os
import shutil
from tempfile import TemporaryDirectory
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore', module='stsci.tools')
from stistools.crrej_from_raw import crrej_from_raw, create_new_crr, HIST_LINE


CRDS_REF_FILE = 'test_crr.fits'
ENVIRON_VAR = 'MOCK_CRDS'


def make_mock_datafile(filename, crds_ref):
    '''Writes a mock RAW datasets with two sets of [SCI, ERR, DQ] exts.
    '''
    hdu = fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(), fits.BinTableHDU(), fits.BinTableHDU(),
        fits.BinTableHDU(), fits.BinTableHDU(), fits.BinTableHDU(), ])

    hdu[0].header['ROOTNAME'] = 'rootname'
    hdu[0].header['DETECTOR'] = 'CCD'
    hdu[0].header['CRREJTAB'] = crds_ref
    hdu[1].header['EXPTIME'] = 1.
    hdu[4].header['EXPTIME'] = 1.

    hdu.writeto(filename, overwrite=True)


def make_mock_crr(filename):
    '''Writes a mock CRREJTAB file with only one row.
    '''
    crr_values = {
        'CRSPLIT':  2,
        'MEANEXP':  4.4,
        'SCALENSE': 0.,
        'INITGUES': 'minimum',
        'SKYSUB':   'mode',
        'CRSIGMAS': '5.0',
        'CRRADIUS': 1.5,
        'CRTHRESH': 0.8,
        'BADINPDQ': 0,
        'CRMASK':   True, }

    new_crr = create_new_crr(crr_values)
    new_crr.write(filename)


class Test_crrej_from_raw:
    @classmethod
    def setup_class(cls):
        cls.refdir = TemporaryDirectory(prefix='test_crrej_from_raw_ref_')
        os.environ[ENVIRON_VAR] = cls.refdir.name + os.path.sep
        make_mock_crr(os.path.join(cls.refdir.name, CRDS_REF_FILE))
        cls.crds_ref_file = f"${ENVIRON_VAR}/{CRDS_REF_FILE}"

    @classmethod
    def teardown_class(cls):
        cls.refdir.cleanup()
        del os.environ[ENVIRON_VAR]

    def setup_method(self, method):
        self.cwd = os.getcwd()
        self.directory = TemporaryDirectory(prefix='test_crrej_from_raw_')
        os.chdir(self.directory.name)
        self.filename = 'test.fits'
        make_mock_datafile(self.filename, crds_ref=self.crds_ref_file)

    def teardown_method(self, method):
        self.directory.cleanup()
        os.chdir(self.cwd)

    def test_user_crrejtab(self):
        '''Checks that the user-specified CRREJTAB overrides the one found in the FITS
        file header.
        '''
        crr = 'local_crr.fits'
        make_mock_crr(crr)
        with fits.open(crr, 'update') as f:
            f[1].data['MEANEXP'][0] = 10000.
        with fits.open(self.filename, 'update') as f:
            f[1].header['EXPTIME'] = 110.  # Out of range in standard test CRRJECTAB
            f[4].header['EXPTIME'] = 110.  # Out of range in standard test CRRJECTAB
        assert crrej_from_raw(self.filename, crrejtab=crr) == 0

    @pytest.mark.parametrize('meanexp,expected_meanexp', [(1.2, 4.4), (4.4, 4.4)])
    def test_meanexp_selection(self, meanexp, expected_meanexp):
        with fits.open(self.filename, 'update') as f:
            f[1].header['EXPTIME'] = meanexp
            f[4].header['EXPTIME'] = meanexp
        assert crrej_from_raw(self.filename) == 0

    def test_bad_meanexp_selection(self):
        with fits.open(self.filename, 'update') as f:
            f[1].header['EXPTIME'] = 110.  # Out of range in test CRRJECTAB
            f[4].header['EXPTIME'] = 110.  # Out of range in test CRRJECTAB
        with pytest.raises(ValueError):
            crrej_from_raw(self.filename)

    def test_scalense(self):
        crrej_from_raw(self.filename, scalense='2.1')
        crrej_from_raw(self.filename, scalense=2.2)

    @pytest.mark.parametrize('initgues', ['minimum', 'median'])
    def test_initgues(self, initgues):
        assert crrej_from_raw(self.filename, initgues=initgues) == 0

    def test_bad_initgues(self):
        with pytest.raises(ValueError):
            crrej_from_raw(self.filename, initgues='bad_value')

    @pytest.mark.parametrize('skysub', ['none', 'mode'])
    def test_skysub(self, skysub):
        assert crrej_from_raw(self.filename, skysub=skysub) == 0

    def test_bad_skysub(self):
        with pytest.raises(ValueError):
            crrej_from_raw(self.filename, skysub='bad_value')

    @pytest.mark.parametrize('crsigmas', ['3', '3,4', '3,4.0', 5, 5.2])
    def test_crsigmas(self, crsigmas):
        crrej_from_raw(self.filename, crsigmas=crsigmas)

    def test_crradius(self):
        '''Checks a user-specified `crradius` value and checks that the
        input FITS file header is returned to its original form.
        '''
        assert crrej_from_raw(self.filename, crradius=2., verbose=True) == 0
        # Check header keywords in input file:
        assert fits.getval(self.filename, ext=0, keyword='CRREJTAB') == \
            f"${ENVIRON_VAR}/{CRDS_REF_FILE}"
        with pytest.raises(KeyError):
            fits.getval(self.filename, ext=0, keyword='PREV_CRR')
        # Check for erroneous presense of HISTORY line in the RAW file:
        try:
            h = list(fits.getval(self.filename, ext=0, keyword='HISTORY'))
        except KeyError:
            h = []
        assert HIST_LINE.strip() not in [x.strip() for x in h]

    @pytest.mark.parametrize('crthresh', [1.1, 3, 0, -2.])
    def test_crthresh(self, crthresh):
        assert crrej_from_raw(self.filename, crthresh=crthresh) == 0

    @pytest.mark.parametrize('badinpdq', [18, 2.0001, '9600'])
    def test_badinpdq(self, badinpdq):
        assert crrej_from_raw(self.filename, badinpdq=badinpdq) == 0

    @pytest.mark.parametrize('crmask', [True, False, None, 1, 0])
    def test_crmask(self, crmask):
        assert crrej_from_raw(self.filename, crmask=crmask) == 0

    def test_input_file_unchanged(self):
        original_file = 'orig_test.fits'
        shutil.copy(self.filename, original_file)
        crrej_from_raw(self.filename, crsigmas='3')
        diff = fits.FITSDiff(original_file, self.filename).report().split('\n')
        assert [x for x in diff if x.strip()][-1] == 'No differences found.', \
            'Input file modified'


class Test_versions:
    def test_version(self):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            crrej_from_raw(None, print_version=True)
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0

    def test_revision(self):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            crrej_from_raw(None, print_revision=True)
        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 0
