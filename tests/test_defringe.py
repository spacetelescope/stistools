from stistools.defringe import normspflat, prepspec, mkfringeflat, defringe
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestDefringe(BaseSTIS):

    input_loc = 'defringe'


    def test_normspflat_g750l(self):
        """Compare normspflat output for a g750l spectrum"""

        sci_file = 'o49x18010'
        flat_file = 'o49x18020'
        output = flat_file+"_nsp_out.fits"

        for input_file in [flat_file+'_raw.fits', sci_file+'_wav.fits']:
            self.get_input_file("input", input_file)
        
        normspflat(flat_file+"_raw.fits",output, 
                   do_cal=True,wavecal=sci_file+"_wav.fits")

        outputs = [(output, flat_file+"_nsp.fits")]
        self.compare_outputs(outputs)


    def test_prepspec_g750l(self):
        """Compare prepspec output for a g750l spectrum"""
        
        sci_file = 'o49x18010'
        flat_file = 'o49x18020'
        output = sci_file+"_crj.fits"

        for input_file in [sci_file+'_raw.fits', sci_file+'_wav.fits']:
            self.get_input_file("input", input_file)

        prepspec(sci_file+"_raw.fits")

        outputs = [(output, sci_file+"_crj.fits")]
        self.compare_outputs(outputs)

        
    def test_mkfringeflat_g750l(self):
        """compare mkfringeflat output for a g750l spectrum"""

        sci_file = 'o49x18010'
        flat_file = 'o49x18020'
        output = flat_file+"_frr_out.fits"

        for input_file in [sci_file+'_crj.fits', flat_file+'_nsp.fits']:
            self.get_input_file("input", input_file)

        mkfringeflat(sci_file+"_crj.fits", flat_file+"_nsp.fits", output)

        outputs = [(output, flat_file+"_frr.fits")]
        self.compare_outputs(outputs)


    def test_defringe_g750l(self):
        """compare defringe output for a g750l spectrum"""
        
        sci_file = 'o49x18010'
        flat_file = 'o49x18020'

        for input_file in [sci_file+'_crj.fits', flat_file+'_frr.fits']:
            self.get_input_file("input", input_file)

        output = defringe(f"{sci_file}_crj.fits",f"{flat_file}_frr.fits")

        outputs = [(output, sci_file+"_drj.fits")]
        self.compare_outputs(outputs)


    def test_normspflat_g750m(self):
        """compare normspflat output for a g750m spectrum"""
        
        sci_file = "oe36m10g0"
        flat_file = "oe36m10j0"
        output = flat_file+"_frr_out.fits"

        for input_file in [flat_file+'_raw.fits', sci_file+'_wav.fits']:
            self.get_input_file("input", input_file)
        
        normspflat(flat_file+"_raw.fits", output,
                   do_cal=True, wavecal=sci_file+"_wav.fits")

        outputs = [(output, flat_file+"_nsp.fits")]
        self.compare_outputs(outputs)


    def test_prepspec_g750m(self):
        """Compare prepspec output for a g750m spectrum"""

        sci_file = "oe36m10g0"
        flat_file = "oe36m10j0"
        output = sci_file+"_sx2.fits"

        for input_file in [sci_file+'_raw.fits', sci_file+'_wav.fits']:
            self.get_input_file("input", input_file)

        prepspec(sci_file+"_raw.fits")

        outputs = [(output, sci_file+"_sx2.fits")]
        self.compare_outputs(outputs)


    def test_mkfringeflat_g750m(self):
        """compare mkfringeflat output for a g750m spectrum"""
        
        sci_file = "oe36m10g0"
        flat_file = "oe36m10j0"

        output = flat_file+"_frr_out.fits"

        for input_file in [sci_file+'_sx2.fits', flat_file+'_nsp.fits']:
            self.get_input_file("input", input_file)

        mkfringeflat(sci_file+"_sx2.fits", flat_file+"_nsp.fits", output,
                     beg_shift=-1.0, end_shift=0.5, shift_step=0.1,
                     beg_scale=0.8, end_scale=1.5, scale_step=0.04)

        outputs = [(output, flat_file+"_frr.fits")]
        self.compare_outputs(outputs)


    def test_defringe_g750m(self):
        """compare defringe output for a g750m spectrum"""

        sci_file = "oe36m10g0"
        flat_file = "oe36m10j0"

        for input_file in [sci_file+'_sx2.fits', flat_file+'_frr.fits']:
            self.get_input_file("input", input_file)

        output = defringe(f"{sci_file}_sx2.fits", f"{flat_file}_frr.fits")

        outputs = [(output, sci_file+"_s2d.fits")]
        self.compare_outputs(outputs)
