from stistools.poisson_err import poisson_err
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestInttag(BaseSTIS):

    input_loc = 'poisson_err'

    def test_accum_lores(self):
        """Compare output for a single-order FUV spectrum"""
        self.get_data("input", "obgh07020_x1d.fits")
        output = "obgh07020_PCI_x1d_out.fits"
        poisson_err("obgh07020_x1d.fits", output)

        outputs = [(output, "obgh07020_PCI_x1d.fits")]
        self.compare_outputs(outputs)

    def test_accum_hires(self):
        """Compare output for an NUV echelle mode spectrum"""
        self.get_data("input", "oep502040_x1d.fits")
        output = "oep502040_PCI_x1d_out.fits"
        poisson_err("oep502040_x1d.fits", output)

        outputs = [(output, "oep502040_PCI_x1d.fits")]
        self.compare_outputs(outputs)
