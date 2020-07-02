import pytest

from stistools.wx2d import wx2d
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestWx2d(BaseSTIS):

    input_loc = 'wx2d'
    ref_loc = 'wx2d'


    def test_wx2d_t1(self):
        """
        Test for  wx2d using rows parameter
        """

        # Prepare input files.
        self.get_input_file("../stisnoise/input", "o6ih10060_crj.fits")

        # Run wx2d
        wx2d("o6ih10060_crj.fits", "o6ih10060_wx2d.fits",
             wavelengths="o6ih10060_wl.fits", helcorr="perform",
             rows=(843, 942))

        # Compare results
        outputs = [("o6ih10060_wx2d.fits", "o6ih10060_wx2d_ref.fits"),
                   ("o6ih10060_wl.fits", "o6ih10060_wl_ref.fits")]
        self.compare_outputs(outputs, delete_history=True)


    def test_wx2d_t2(self):
        """
        Test for wx2d
        """

        # Prepare input files.
        self.get_input_file("input", "o4d301030_crj.fits")

        # Run wx2d
        wx2d("o4d301030_crj.fits", "o4d301030_wx2d.fits", helcorr="perform")

        # Compare results
        # Compare results
        outputs = [("o4d301030_wx2d.fits", "o4d301030_wx2d_ref.fits")]
        self.compare_outputs(outputs, delete_history=True)
