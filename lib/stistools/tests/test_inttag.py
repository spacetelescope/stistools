from ..inttag import inttag
from astropy.io import fits
import numpy as np
import pytest

threshold = 0.0001  # Acceptable ratio of different pixels to identical pixels

def test_accum_lores():
    """Compare accum image output for a single lowres imset"""
    iraf_file = 'data/inttag/lores.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/oddv01050_tag.fits", output, highres=False)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data

    assert (np.sum(output_data != iraf_data) / np.sum(output_data == iraf_data) < threshold) == True


def test_accum_hires():
    """Compare accum image output for a single highres imset"""
    iraf_file = 'data/inttag/hires.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/oddv01050_tag.fits", output, highres=True)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data

    assert (np.sum(output_data != iraf_data) / np.sum(output_data == iraf_data) < threshold) == True


def test_accum_gtigap():
    """Compare accum image output for a single imset with a GTI gap"""
    iraf_file = 'data/inttag/gtigap.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/gtigap_tag.fits", output, highres=False)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data

    assert (np.sum(output_data != iraf_data) / np.sum(output_data == iraf_data) < threshold) == True


def test_accum_allevents():
    """Compare accum image output in allevents mode"""
    iraf_file = 'data/inttag/allevents.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/gtigap_tag.fits", output, highres=False, allevents=True)
    iraf_data = fits.open(iraf_file)[1].data
    output_data = fits.open(output)[1].data

    assert (np.sum(output_data != iraf_data) / np.sum(output_data == iraf_data) < threshold) == True


def test_accum_rcount():
    """Compare accum image multi-imset output"""
    iraf_file = 'data/inttag/rcount.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/gtigap_tag.fits", output, starttime=700, increment=200, rcount=5, highres=False)
    n_sci = int(fits.open("data/inttag/rcount.fits")[0].header['NEXTEND'] / 3)
    for sci in range(n_sci):
        iraf_data = fits.open(iraf_file)['SCI', sci+1].data
        output_data = fits.open(output)['SCI', sci+1].data

        assert (np.sum(output_data != iraf_data) / np.sum(output_data == iraf_data) < threshold) == True