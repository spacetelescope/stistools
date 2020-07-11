from stistools.inttag import inttag
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestInttag(BaseSTIS):

    input_loc = 'inttag'

    def test_accum_lores(self):
        """Compare accum image output for a single lowres imset"""
        self.get_data("input", "oddv01050_tag.fits")
        output = "inttag_accum_lores_out.fits"
        inttag("oddv01050_tag.fits", output, highres=False)

        outputs = [(output, "inttag_accum_lores.fits")]
        self.compare_outputs(outputs)

    def test_accum_hires(self):
        """Compare accum image output for a single highres imset"""
        self.get_data("input", "oddv01050_tag.fits")
        output = "inttag_accum_hires_out.fits"
        inttag("oddv01050_tag.fits", output, highres=True)

        outputs = [(output, "inttag_accum_hires.fits")]
        self.compare_outputs(outputs)

    def test_accum_gtigap(self):
        """Compare accum image output for a single imset with a GTI gap"""
        self.get_data("input", "gtigap_tag.fits")
        output = "inttag_accum_gtigap_out.fits"
        inttag("gtigap_tag.fits", output, highres=False)

        outputs = [(output, "inttag_accum_gtigap.fits")]
        self.compare_outputs(outputs)

    def test_accum_allevents(self):
        """Compare accum image output in allevents mode"""
        self.get_data("input", "gtigap_tag.fits")
        output = "inttag_accum_allevents_out.fits"
        inttag("gtigap_tag.fits", output, highres=False, allevents=True)

        outputs = [(output, "inttag_accum_allevents.fits")]
        self.compare_outputs(outputs)

    def test_accum_rcount(self):
        """Compare accum image multi-imset output"""
        self.get_data("input", "gtigap_tag.fits")
        output = "inttag_accum_rcount_out.fits"
        inttag("gtigap_tag.fits", output, starttime=700, increment=200, rcount=5, highres=False)

        outputs = [(output, "inttag_accum_rcount.fits")]
        self.compare_outputs(outputs)

    def test_primary_hdr_nogap(self):
        """Compare generated primary header keywords (with no gti gap in data)"""
        self.get_data("input", "ob3001xqq_tag.fits")
        output = "inttag_prihdr_nogap_out.fits"
        inttag("ob3001xqq_tag.fits", output)

        outputs = [(output, "inttag_prihdr_nogap.fits")]
        self.compare_outputs(outputs)

    def test_primary_hdr_gap(self):
        """Compare generated primary header keywords (with some gti gap in data)"""
        self.get_data("input", "gtigap_tag.fits")
        output = "inttag_prihdr_gap_out.fits"
        inttag("gtigap_tag.fits", output)

        outputs = [(output, "inttag_prihdr_gap.fits")]
        self.compare_outputs(outputs)

    def test_exptime_truncation(self):
        """Check if inttag handles input rcount sizes correctly by truncating at the last event"""
        self.get_data("input", "od7s08010_tag.fits")
        output = "inttag_exptime_trunc_out.fits"
        inttag("od7s08010_tag.fits", output, increment=1675., rcount=2, starttime=0., allevents=False)

        outputs = [(output, "inttag_exptime_trunc.fits")]
        self.compare_outputs(outputs)

