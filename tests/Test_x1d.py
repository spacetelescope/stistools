from stistools.x1d import x1d
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestX1d(BaseSTIS):

    input_loc = 'x1d'
    ref_loc = 'x1d/ref'

    def test_x1d(self):
        """
        Basic x1d test, mostly using default parameters
        """

        # Prepare input files.
        self.get_input_file("input", "o56j02020_flt.fits")

        # Run basic2d
        x1d("o56j02020_flt.fits", output="x1d_lev3.fits", ctecorr="omit")

        # Compare results
        outputs = [("x1d_lev3.fits", "x1d_lev3_ref.fits")]
        self.compare_outputs(outputs)
