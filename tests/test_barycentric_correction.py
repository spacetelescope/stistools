#!/usr/bin/env python3
from astropy.io import fits
from stistools.barycentric_correction import barycentric_correction, OrbFileError
from .resources import BaseSTIS
import pytest

# TO DO:  Download one local copy of the ORB_FILE for the entire suite.
#         Be careful about modifications due to `test_orbile_column_case` and
#         `test_orbfile_column_missing`.
#         This does not currently work because class methods can't find the
#         proper data location on Artifactory (`BaseSTIS.get_data()`).


ORB_FILE = 'p2o0000r.fit'  # 2017-02-24 00:00:00 - 2017-02-24 00:00:00

STIS_FILES = [
    'od9m97020_raw.fits',
    'od9m97020_tag.fits',
    'od9m97020_x1d.fits']

COS_FILES = [
    'ldad10m0q_rawtag_a.fits',
    'ldad10m0q_corrtag_a.fits',
    'ldad10m0q_x1d.fits',
    'ldad10020_x1dsum.fits']


@pytest.mark.bigdata
@pytest.mark.slow
class TestBarycentricCorrection(BaseSTIS):

    input_loc = 'barycentric_correction'

    @pytest.mark.parametrize('filename', STIS_FILES)
    def test_stis_jpl(self, filename):
        '''Compare output for STIS files using JPL ephemeris.
        '''
        self.get_data('input', filename)
        barycentric_correction(filename)
        self.compare_outputs([(filename, 'jpl_' + filename)])

    @pytest.mark.parametrize('filename', STIS_FILES)
    def test_stis_orbfile(self, filename):
        '''Compare output for STIS files using HST orbfile ephemeris.
        '''
        self.get_data('input', filename)
        self.get_data('input', ORB_FILE)
        barycentric_correction(filename, hst_orb=ORB_FILE)
        self.compare_outputs([(filename, 'hst_orb_' + filename)])

    @pytest.mark.parametrize('filename', COS_FILES)
    def test_cos_jpl(self, filename):
        '''Compare output for COS files using JPL ephemeris.
        '''
        self.get_data('input', filename)
        barycentric_correction(filename)
        self.compare_outputs([(filename, 'jpl_' + filename)])

    @pytest.mark.parametrize('filename', COS_FILES)
    def test_cos_orbfile(self, filename):
        '''Compare output for COS files using HST orbfile ephemeris.
        '''
        self.get_data('input', filename)
        self.get_data('input', ORB_FILE)
        barycentric_correction(filename, hst_orb=ORB_FILE)
        self.compare_outputs([(filename, 'hst_orb_' + filename)])

    @pytest.mark.parametrize('filename', STIS_FILES[1:])  # tag, x1d
    def test_stis_utc(self, filename):
        '''Test calculation of UTC time system values.
        '''
        self.get_data('input', filename)
        barycentric_correction(filename, time_system='UTC', time_script=True)
        self.compare_outputs([(filename, 'jpl_utc_' + filename)])

    def test_orbile_column_case(self):
        '''Rename orb_file's column from 'Time' to 'TIME' to show case-insensitivity.
        '''
        test_file = STIS_FILES[0]  # raw
        self.get_data('input', test_file)
        self.get_data('input', ORB_FILE)

        # Modify orb_file to rename column:
        #   'p2o0000r.fit' has a 'Time' column.  Try 'TIME'...
        with fits.open(ORB_FILE, 'update') as f:
            f[1].data.columns.change_name('Time', 'TIME')

        barycentric_correction(test_file, hst_orb=ORB_FILE)
        self.compare_outputs([(test_file, 'hst_orb_' + test_file)])

    def test_orbfile_column_missing(self):
        '''Rename orb_file's column from 'Time' to 'BAD' to enforce ValueError is raised.
        '''
        test_file = STIS_FILES[0]  # raw
        self.get_data('input', test_file)
        self.get_data('input', ORB_FILE)

        with fits.open(ORB_FILE, 'update') as f:
            f[1].data.columns.change_name('Time', 'BAD')

        with pytest.raises(ValueError):
            barycentric_correction(test_file, hst_orb=ORB_FILE)

    def test_orbfile_exception(self):
        '''Test running the correction outside of the hst_orb file's date range.
        '''
        test_file = 'oda913010_x1d.fits'  # start time:  2017-02-27 01:06:28
        self.get_data('input', test_file)
        self.get_data('input', ORB_FILE)

        with pytest.raises(OrbFileError):
            barycentric_correction(test_file, hst_orb=ORB_FILE)
