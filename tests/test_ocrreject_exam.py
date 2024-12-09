from stistools.ocrreject_exam import ocrreject_exam
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestOcrrejectExam(BaseSTIS):

    input_loc = 'ocrreject_exam'
    ref_loc = 'ocrreject_exam/ref'

    input_list = ["odvkl1040_flt.fits", "odvkl1040_sx1.fits"]

    # Make input file string
    input_file_string = ", ".join(input_list)

    def test_ocrrject_exam(self):
        """
        This regression test runs the task on the known problematic dataset odvkl1040.
        The resulting output dictionary  is compared to a
        reference file using 'FITSDIFF'. COULD I JUST PULL OUT THE KEY/VAL PAIRS AND CHECK IF THEY'RE SIMILAR? 
        """

        # Prepare input files.
        for filename in self.input_list:
            self.get_input_file("input", filename)

        # Run ocrreject_exam
        ocrreject(self.input_file_string, output="ocrreject_lev3_crj.fits")

        # Compare results
        outputs = [("ocrreject_lev3_crj.fits", "ocrreject_lev3_crj_ref.fits")]
        self.compare_outputs(outputs)