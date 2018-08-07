from ..calstis import calstis
from .resources import BaseSTIS
from .helpers.mark import require_bigdata

from astropy.io import fits

@require_bigdata
class TestSTIS(BaseSTIS):

    def test_ccd_imaging(self):
        """ This test is for calstis on CCD imaging data
        """

        # Prepare input files.
        input_file = self.get_input_file("input", "oci402010_raw.fits")
        input_epc1 = self.get_data("input", "oci402dnj_epc.fits")
        input_epc2 = self.get_data("input", "oci402doj_epc.fits")

        outroot = "stis_test1"
        outfile = outroot + "_flt.fits"
        outcrj = outroot + "_crj.fits"
        outsx2 = outroot + "_sx2.fits"
        reffile_flt = 'reference_stis_01_flt.fits'
        reffile_crj = 'reference_stis_01_crj.fits'
        reffile_sx2 = 'reference_stis_01_sx2.fits'

        calstis("oci402010_raw.fits", outroot=outroot)

        #fail dang you
        fits.setval('stis_test1_flt.fits', 'DARKCORR', value="NA")

        # Compare results
        outputs = [(outfile, reffile_flt), (outcrj, reffile_crj),
                   (outsx2, reffile_sx2)]
        self.compare_outputs(outputs)
