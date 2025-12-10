from stistools.mktrace import mktrace
from .resources import BaseSTIS
import re
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestMktrace(BaseSTIS):

    input_loc = 'mktrace'
    ref_loc = 'mktrace'

    atol = 1e-14
    rtol = 2e-4

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
        m = re.search(
            r"Traces were rotated by ([0-9.]+) degrees.*?"
            r"trace is centered on row ([0-9.]+)",
            captured.out,
            re.S)
        angle = float(m.group(1))
        row = float(m.group(2))

        assert angle == pytest.approx(0.0066273851, rel=1e-6)
        assert row == pytest.approx(509.2475494292, rel=1e-6))

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
        m = re.search(
            r"Traces were rotated by ([0-9.]+) degrees.*?"
            r"trace is centered on row ([0-9.]+)",
            captured.out,
            re.S)
        angle = float(m.group(1))
        row = float(m.group(2))

        assert angle == pytest.approx(0.0317182356, rel=1e-6)
        assert row == pytest.approx(893.7059707738, rel=1e-6)

        outputs = [("o8pp31020_crj_1dt.fits", "o8pp31020_crj_1dt_ref.fits"),
                   ("o8pp31020_crj_1dt_interp.fits", "o8pp31020_crj_1dt_interp_ref.fits"),
                   ("o8pp31020_crj_1dt_interpfit.fits", "o8pp31020_crj_1dt_interpfit_ref.fits"),
                   ("o8pp31020_crj_1dt_sci.fits", "o8pp31020_crj_1dt_sci_ref.fits"),
                   ("o8pp31020_crj_1dt_scifit.fits", "o8pp31020_crj_1dt_scifit_ref.fits")]
        self.compare_outputs(outputs)
