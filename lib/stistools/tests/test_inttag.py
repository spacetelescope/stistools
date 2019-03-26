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


def test_primary_hdr_nogap():
    """Compare generated primary header keywords (with no gti gap in data)"""
    iraf_file = 'data/inttag/ob3001xqq_raw.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/ob3001xqq_tag.fits", output)
    iraf_hdr = fits.open(iraf_file)[0].header
    output_hdr = fits.open(output)[0].header

    assert output_hdr['NEXTEND'] == iraf_hdr['NEXTEND']
    assert output_hdr['TDATEOBS'] == iraf_hdr['TDATEOBS']
    assert output_hdr['TTIMEOBS'] == iraf_hdr['TTIMEOBS']
    assert output_hdr['TEXPTIME'] == pytest.approx(iraf_hdr['TEXPTIME'], 0.1)
    assert output_hdr['TEXPSTRT'] == pytest.approx(iraf_hdr['TEXPSTRT'], 0.1)
    assert output_hdr['TEXPEND'] == pytest.approx(iraf_hdr['TEXPEND'], 0.1)
    assert output_hdr['LORSCORR'] == "COMPLETE"


def test_primary_hdr_gap():
    """Compare generated primary header keywords (with some gti gap in data)"""
    iraf_file = 'data/inttag/gtigap.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/gtigap_tag.fits", output)
    iraf_hdr = fits.open(iraf_file)[0].header
    output_hdr = fits.open(output)[0].header

    assert output_hdr['NEXTEND'] == iraf_hdr['NEXTEND']
    assert output_hdr['TDATEOBS'] == iraf_hdr['TDATEOBS']
    assert output_hdr['TTIMEOBS'] == iraf_hdr['TTIMEOBS']
    assert output_hdr['TEXPTIME'] == pytest.approx(iraf_hdr['TEXPTIME'], 0.1)
    assert output_hdr['TEXPSTRT'] == pytest.approx(iraf_hdr['TEXPSTRT'], 0.1)
    assert output_hdr['TEXPEND'] == pytest.approx(iraf_hdr['TEXPEND'], 0.1)


def test_science_hdr_lores():
    """Compare generated science header keywords in lores mode"""
    iraf_file = 'data/inttag/lores.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/oddv01050_tag.fits", output)
    iraf_hdr = fits.open(iraf_file)[1].header
    output_hdr = fits.open(output)[1].header

    assert output_hdr['WCSAXES'] == iraf_hdr['WCSAXES']
    assert output_hdr['LTM1_1'] == iraf_hdr['LTM1_1']
    assert output_hdr['LTM2_2'] == iraf_hdr['LTM2_2']
    assert output_hdr['LTV1'] == iraf_hdr['LTV1']
    assert output_hdr['LTV2'] == iraf_hdr['LTV2']
    assert output_hdr['CD1_1'] == iraf_hdr['CD1_1']
    assert output_hdr['CD1_2'] == iraf_hdr['CD1_2']
    assert output_hdr['CD2_1'] == iraf_hdr['CD2_1']
    assert output_hdr['CD2_2'] == iraf_hdr['CD2_2']
    assert output_hdr['CRPIX1'] == iraf_hdr['CRPIX1']
    assert output_hdr['CRPIX2'] == iraf_hdr['CRPIX2']


def test_science_hdr_hires():
    """Compare generated science header keywords in lores mode"""
    iraf_file = 'data/inttag/hires.fits'
    output = "data/inttag/test.fits"
    inttag("data/inttag/oddv01050_tag.fits", output, highres=True)
    iraf_hdr = fits.open(iraf_file)[1].header
    output_hdr = fits.open(output)[1].header

    assert output_hdr['WCSAXES'] == iraf_hdr['WCSAXES']
    assert output_hdr['LTM1_1'] == iraf_hdr['LTM1_1']
    assert output_hdr['LTM2_2'] == iraf_hdr['LTM2_2']
    assert output_hdr['LTV1'] == iraf_hdr['LTV1']
    assert output_hdr['LTV2'] == iraf_hdr['LTV2']
    assert output_hdr['CD1_1'] == iraf_hdr['CD1_1']
    assert output_hdr['CD1_2'] == iraf_hdr['CD1_2']
    assert output_hdr['CD2_1'] == iraf_hdr['CD2_1']
    assert output_hdr['CD2_2'] == iraf_hdr['CD2_2']
    assert output_hdr['CRPIX1'] == iraf_hdr['CRPIX1']
    assert output_hdr['CRPIX2'] == iraf_hdr['CRPIX2']
