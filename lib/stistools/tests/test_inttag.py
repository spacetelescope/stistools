from ..inttag import inttag
from astropy.io import fits
import numpy as np


def iraf_compare_accum_lores():
    """Compare accum image output for a single lowres imset"""
    iraf_file = 'data/lores.fits'
    output = "test.fits"
    inttag("data/oddv01050_tag.fits", output, highres=False)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data


def iraf_compare_accum_hires():
    """Compare accum image output for a single highres imset"""
    iraf_file = 'data/hires.fits'
    output = "test.fits"
    inttag("data/oddv01050_tag.fits", output, highres=True)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data


def iraf_compare_accum_gtigap():
    """Compare accum image output for a single imset with a GTI gap"""
    iraf_file = 'data/gtigap.fits'
    output = "test.fits"
    inttag("data/gtigap_tag.fits", output, highres=False)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data


def iraf_compare_accum_allevents():
    """Compare accum image output in allevents mode"""
    iraf_file = 'data/allevents.fits'
    output = "test.fits"
    inttag("data/gtigap_tag.fits", output, highres=False, allevents=True)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data


def iraf_compare_accum_rcount():
    """Compare accum image multi-imset output"""
    iraf_file = 'data/rcount.fits'
    output = "test.fits"
    inttag("data/gtigap_tag.fits", output, starttime=700, increment=200, rcount=5, highres=False)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data