#!/usr/bin/env python3
from stistools.barycentric_correction import bary_corr
from .resources import BaseSTIS
import pytest


@pytest.mark.bigdata
@pytest.mark.slow
class TestInttag(BaseSTIS):

    input_loc = 'barycentric_correction'

    def test_raw_JPL(self):
        """Compare output for a x1d file."""
        self.get_data("input", "od9m97020_raw.fits".fits")
        #output = "od9m97020_barycorr_raw_out.fits"
        
        
        bary_corr("od9m97020_raw.fits")

        outputs = [(output, "od9m97020_raw_JPL_out.fits")]
        self.compare_outputs(outputs)
        
    def test_raw_orbfile(self):
        """Compare output for a x1d file."""
        self.get_data("input", "od9m97020_raw.fits".fits")
        #output = "od9m97020_barycorr_raw_out.fits"
        
        
        bary_corr("od9m97020_raw.fits", hst_orb="p2o0000r_TIME.fit")

        outputs = [(output, "od9m97020_raw_orbfile_out.fits")]
        self.compare_outputs(outputs)