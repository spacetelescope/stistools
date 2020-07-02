from stistools.calstis import calstis
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestCalstis(BaseSTIS):

    input_loc = 'calstis'
    ref_loc = 'calstis/ref'

    def test_ccd_imaging(self):
        """
        This test is for calstis on CCD imaging data
        """

        # Prepare input files.
        self.get_input_file("input", "oci402010_raw.fits")
        self.get_data("input", "oci402dnj_epc.fits")
        self.get_data("input", "oci402doj_epc.fits")

        outroot = "stis_test1"
        outfile = outroot + "_flt.fits"
        outcrj = outroot + "_crj.fits"
        outsx2 = outroot + "_sx2.fits"
        reffile_flt = 'reference_stis_01_flt.fits'
        reffile_crj = 'reference_stis_01_crj.fits'
        reffile_sx2 = 'reference_stis_01_sx2.fits'

        calstis("oci402010_raw.fits", outroot=outroot)

        # Compare results
        outputs = [(outfile, reffile_flt), (outcrj, reffile_crj),
                   (outsx2, reffile_sx2)]
        self.compare_outputs(outputs)

    # def test1_lev1_FUV(self):
    #     """
    #     This test is for level 1? FUV data
    #     """
    #
    #     # Prepare input files.
    #     self.get_input_file("input", "o5cl02040_raw.fits")
    #     self.get_input_file("input", "o5cl02040_wav.fits")
    #     outroot = "calstis_lev1_FUVspec"
    #
    #     # Run test
    #     calstis("o5cl02040_raw.fits", outroot=outroot)
    #
    #     # Compare results
    #     outputs = [(outroot+"_flt.fits", outroot+"_flt_ref.fits"),
    #                (outroot+"_x1d.fits", outroot+"_x1d_ref.fits")]
    #     self.compare_outputs(outputs)


    def test2_lev1_FUV(self):
        """
        This test is for level 1? FUV data
        """

        # Prepare input files.
        self.get_input_file("input", "odj102010_raw.fits")
        self.get_input_file("input", "odj102010_wav.fits")
        outroot = "calstis_1lev1_FUVspec"

        # Run test
        calstis("odj102010_raw.fits", outroot=outroot)

        # Compare results
        outputs = [(outroot+"_flt.fits", outroot+"_flt_ref.fits"),
                   (outroot+"_x1d.fits", outroot+"_x1d_ref.fits")]
        self.compare_outputs(outputs)

    def test_lev2_CCD(self):
        """
        This test is for level 2 CCD data
        """

        # Prepare input files.
        self.get_input_file("input", "o3wd01060_raw.fits")
        self.get_input_file("input", "o3wd01060_wav.fits")
        self.get_data("input", "o3wd01060_spt.fits")
        outroot = "calstis_lev2_CCD"

        # Run test
        calstis("o3wd01060_raw.fits", outroot=outroot)

        # Compare results
        outputs = [(outroot + "_flt.fits", outroot + "_flt_ref.fits"),
                   (outroot + "_crj.fits", outroot + "_crj_ref.fits"),
                   (outroot + "_sx1.fits", outroot + "_sx1_ref.fits"),
                   (outroot + "_sx2.fits", outroot + "_sx2_ref.fits")]
        self.compare_outputs(outputs)

    def test_lev3_FUV(self):
        """
        This test is for level 3 FUV data
        """

        # Prepare input files.
        self.get_input_file("input", "o5in01tnq_raw.fits")
        outroot = "calstis_lev3_FUV"

        # Run test
        calstis("o5in01tnq_raw.fits", outroot=outroot)

        # Compare results
        outputs = [(outroot + "_flt.fits", outroot + "_flt_ref.fits"),
                   (outroot + "_x2d.fits", outroot + "_x2d_ref.fits")]
        self.compare_outputs(outputs)

    def test_lev3_NUV(self):
        """
        This test is for level 3 NUV data
        """

        # Prepare input files.
        self.get_input_file("input", "o6d806030_raw.fits")
        self.get_data("input", "o6d806030_wav.fits")
        outroot = "calstis_lev3_NUV"

        # Run test
        calstis("o6d806030_raw.fits", outroot=outroot)

        # Compare results
        outputs = [(outroot + "_flt.fits", outroot + "_flt_ref.fits"),
                   (outroot + "_x1d.fits", outroot + "_x1d_ref.fits"),
                   (outroot + "_x2d.fits", outroot + "_x2d_ref.fits")]
        self.compare_outputs(outputs)
