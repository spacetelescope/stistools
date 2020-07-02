from stistools.mktrace import mktrace
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestMktrace(BaseSTIS):

    input_loc = 'mktrace'
    ref_loc = 'mktrace'

    atol = 1e-14

    def test_mktrace_t1(self, capsys):
        """
        This tests a basic usage of stis.mktrace. To plot the results run
        plot*.py
        """

        # Prepare input files.
        self.get_input_file("input", "o45a03010_crj.fits")

        capsys.readouterr()

        # Run mktrace
        mktrace('o45a03010_crj.fits')

        # Compare results
        captured = capsys.readouterr()
        assert captured.out == "Traces were rotated by 0.0066273929 degrees \n" \
                               "\n" \
                               "trace is centered on row 509.2475494293\n"

        outputs = [("o45a03010_crj_1dt.fits", "o45a03010_crj_1dt_ref.fits"),
                   ("o45a03010_crj_1dt_interp.fits", "o45a03010_crj_1dt_interp_ref.fits"),
                   ("o45a03010_crj_1dt_interpfit.fits", "o45a03010_crj_1dt_interpfit_ref.fits"),
                   ("o45a03010_crj_1dt_sci.fits", "o45a03010_crj_1dt_sci_ref.fits"),
                   ("o45a03010_crj_1dt_scifit.fits", "o45a03010_crj_1dt_scifit_ref.fits")]
        self.compare_outputs(outputs)

    def test_mktrace_t2(self, capsys):
        """
        This test uses E1 aperture position and places the
        target near row 900. Originally there was a problem
        with the interpolated trace due to incorrect determination
        of the index of the trace inthe trace table.
        """

        # Prepare input files.
        self.get_input_file("input", "o8pp31020_crj.fits")

        capsys.readouterr()

        # Run mktrace
        mktrace('o8pp31020_crj.fits')

        # Compare results
        captured = capsys.readouterr()
        assert captured.out == "Traces were rotated by 0.0317233207 degrees \n" \
                               "\n" \
                               "trace is centered on row 893.7059688464\n"

        outputs = [("o8pp31020_crj_1dt.fits", "o8pp31020_crj_1dt_ref.fits"),
                   ("o8pp31020_crj_1dt_interp.fits", "o8pp31020_crj_1dt_interp_ref.fits"),
                   ("o8pp31020_crj_1dt_interpfit.fits", "o8pp31020_crj_1dt_interpfit_ref.fits"),
                   ("o8pp31020_crj_1dt_sci.fits", "o8pp31020_crj_1dt_sci_ref.fits"),
                   ("o8pp31020_crj_1dt_scifit.fits", "o8pp31020_crj_1dt_scifit_ref.fits")]
        self.compare_outputs(outputs)
