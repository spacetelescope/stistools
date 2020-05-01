from stistools.stisnoise import stisnoise
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestStisnoise(BaseSTIS):

    input_loc = 'stisnoise'
    ref_loc = 'stisnoise'

    # Really should add some output and return array checks for these tests

    def test_stisnoise_t1(self, capsys):
        """
        stisnoise with boxcar test and check output
        """

        # Prepare input files.
        self.get_input_file("input", "o6ih10060_crj.fits")

        capsys.readouterr()

        # Run stisnoise
        stisnoise('o6ih10060_crj.fits', outfile='o6ih10060_fcrj1.fits',
                  boxcar=5)

        # Compare results
        captured = capsys.readouterr()
        assert captured.out == "Target: V1016-CYG, Amp: D, Gain: 1\n"

        outputs = [("o6ih10060_fcrj1.fits", "o6ih10060_fcrj1_ref.fits")]
        self.compare_outputs(outputs)

    def test_stisnoise_t2(self):
        """
        stisnoise with window test
        """

        # Prepare input files.
        self.get_input_file("input", "o6ih10060_crj.fits")

        # Run stisnoise
        stisnoise('o6ih10060_crj.fits', outfile='o6ih10060_fcrj2.fits',
                  window=[0.5, 0.5, 0.1])

        # Compare results
        outputs = [("o6ih10060_fcrj2.fits", "o6ih10060_fcrj2_ref.fits")]
        self.compare_outputs(outputs)
