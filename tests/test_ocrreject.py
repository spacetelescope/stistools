from stistools.ocrreject import ocrreject
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestOcrreject(BaseSTIS):

    input_loc = 'ocrreject'
    ref_loc = 'ocrreject/ref'

    input_list = ["o58i01q7q_flt.fits", "o58i01q8q_flt.fits",
                  "o58i01q9q_flt.fits", "o58i01qcq_flt.fits",
                  "o58i01qdq_flt.fits", "o58i01qeq_flt.fits",
                  "o58i01qhq_flt.fits", "o58i01qiq_flt.fits",
                  "o58i01qjq_flt.fits", "o58i01qmq_flt.fits",
                  "o58i01qnq_flt.fits", "o58i01qoq_flt.fits",
                  "o58i01qrq_flt.fits", "o58i01qsq_flt.fits"]

    # Make input file string
    input_file_string = ", ".join(input_list)

    def test_ocrrject_lev2(self):
        """
        This regression test for this level of this task is different than
        level three in two ways.  Two parameters are set to give an initial
        guess to the sky value and define a sky subraction method. It also
        removes cosmic rays from 14 STIS/CCD images and creates a single
        'clean' image which is compared to a reference file using 'FITSDIFF'.
        """

        # Prepare input files.
        for filename in self.input_list:
            self.get_input_file("input", filename)

        # Run ocrreject
        ocrreject(self.input_file_string, output="ocrreject_lev2_crj.fits",
                  initgues="med", skysub="mode")

        # Compare results
        outputs = [("ocrreject_lev2_crj.fits", "ocrreject_lev2_crj_ref.fits")]
        self.compare_outputs(outputs)

    def test_ocrrject_lev3(self):
        """
        This regression test for this level on this task is a simple default
        parameter execution of the task.  It attempts to remove cosmic rays
        from 14 STIS/CCD images.  The resulting calibration is compared to a
        reference file using 'FITSDIFF'.
        """

        # Prepare input files.
        for filename in self.input_list:
            self.get_input_file("input", filename)

        # Run ocrreject
        ocrreject(self.input_file_string, output="ocrreject_lev3_crj.fits")

        # Compare results
        outputs = [("ocrreject_lev3_crj.fits", "ocrreject_lev3_crj_ref.fits")]
        self.compare_outputs(outputs)
