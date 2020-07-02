from stistools.basic2d import basic2d
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestBasic2d(BaseSTIS):

    input_loc = 'basic2d'
    ref_loc = 'basic2d/ref'

    def test_basic2d_lev1a(self):
        """
        BASIC2D - stsdas/hst_calib/stis/basic2d: level 1a
        This regression test for this level of this task is a running of the
        task with ALL parameters blazing (equal to 'perform').  This includes
        the additional creation of an output test file for bias levels.
        The input data for this test is a STIS NUV_MAMA CCD image.
        """

        # Prepare input files.
        self.get_input_file("input", "o6d806030_raw.fits")

        # Run basic2d
        basic2d("o6d806030_raw.fits", output="basic2d_lev1a_flt.fits")

        # Compare results
        outputs = [("basic2d_lev1a_flt.fits", "basic2d_lev1a_flt_ref.fits")]
        self.compare_outputs(outputs)

    def test_basic2d_lev3a_blev(self):
        """
        The regression test for this level of this task is a running of the
        task with ONLY the bias level parameter parameter set (equal to
        'perform').  This includes the additional creation of an output test
        file for bias levels. The input data for this test is a STIS CCD image.
        """

        # Prepare input files.
        self.get_input_file("input", "o3wd01060_raw.fits")
        self.get_data("input", "o3wd01060_wav.fits")
        self.get_data("input", "o3wd01060_spt.fits")

        # Run basic2d
        basic2d('o3wd01060_raw.fits', output="basic2d_lev3a_blev.fits",
                outblev="basic2d_lev3a_blev.dat", dqicorr="omit",
                doppcorr="omit", lorscorr="omit", glincorr="omit",
                lflgcorr="omit", biascorr="omit", darkcorr="omit",
                flatcorr="omit", photcorr="omit")

        # Compare results
        outputs = [("basic2d_lev3a_blev.fits", "basic2d_lev3a_blev_ref.fits"),
                   ("basic2d_lev3a_blev.dat", "basic2d_lev3a_blev_ref.dat")]
        self.compare_outputs(outputs)

    def test_basic2d_lev3b_blev(self):
        """
        The regression test for this level of this task is a running of the
        task with the bias level parameter and Doppler smoothing correction
        parameters set (equal to 'perform').  This includes the additional
        creation of an output test file for bias levels. The input data for
        this test is a STIS CCD image.
        """

        # Prepare input files.
        self.get_input_file("input", "o3wd01060_raw.fits")
        self.get_data("input", "o3wd01060_wav.fits")
        self.get_data("input", "o3wd01060_spt.fits")

        # Run basic2d
        basic2d('o3wd01060_raw.fits', output="basic2d_lev3b_blev.fits",
                outblev="basic2d_lev3b_blev.dat", dqicorr="omit",
                lorscorr="omit", glincorr="omit", lflgcorr="omit",
                biascorr="omit", darkcorr="omit", flatcorr="omit",
                photcorr="omit")

        # Compare results
        outputs = [("basic2d_lev3b_blev.fits", "basic2d_lev3b_blev_ref.fits"),
                   ("basic2d_lev3b_blev.dat", "basic2d_lev3b_blev_ref.dat")]
        self.compare_outputs(outputs)

    def test_basic2d_lev3c_blev(self):
        """
        The regression test for this level of this task is a running of the
        task with the bias level, Doppler smoothing correction and BIAS image
        sub. parameters set (equal to 'perform').  This includes the additional
        creation of an output test file for bias levels. The input data for
        this test is a STIS CCD image.
        """

        # Prepare input files.
        self.get_input_file("input", "o3wd01060_raw.fits")
        self.get_data("input", "o3wd01060_wav.fits")
        self.get_data("input", "o3wd01060_spt.fits")

        # Run basic2d
        basic2d('o3wd01060_raw.fits', output="basic2d_lev3c_blev.fits",
                outblev="basic2d_lev3c_blev.dat", dqicorr="omit",
                lorscorr="omit", glincorr="omit", lflgcorr="omit",
                darkcorr="omit", flatcorr="omit", photcorr="omit")

        # Compare results
        outputs = [("basic2d_lev3c_blev.fits", "basic2d_lev3c_blev_ref.fits"),
                   ("basic2d_lev3c_blev.dat", "basic2d_lev3c_blev_ref.dat")]
        self.compare_outputs(outputs)

    def test_basic2d_lev3d_blev(self):
        """
        The regression test for this level of this task is a running of the
        task with the bias level, Doppler smoothing correction, BIAS image
        subtraction, and dark subtraction parameters set (equal to 'perform').
        This includes the additional creation of an output test file for bias
        levels. The input data for this test is a STIS CCD image.
        """

        # Prepare input files.
        self.get_input_file("input", "o3wd01060_raw.fits")
        self.get_data("input", "o3wd01060_wav.fits")
        self.get_data("input", "o3wd01060_spt.fits")

        # Run basic2d
        basic2d('o3wd01060_raw.fits', output="basic2d_lev3d_blev.fits",
                outblev="basic2d_lev3d_blev.dat", dqicorr="omit",
                lorscorr="omit", glincorr="omit", lflgcorr="omit",
                flatcorr="omit", photcorr="omit")

        # Compare results
        outputs = [("basic2d_lev3d_blev.fits", "basic2d_lev3d_blev_ref.fits"),
                   ("basic2d_lev3d_blev.dat", "basic2d_lev3d_blev_ref.dat")]
        self.compare_outputs(outputs)
