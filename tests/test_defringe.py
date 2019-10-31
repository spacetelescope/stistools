import os
import tempfile
import numpy as np

from astropy.io import fits

from stistools.defringe import defringe

shape = (5, 7)
datasize = shape[0] * shape[1]

# Name of a temporary directory.  See functions save_file and clean_up.
temp_dir = None
# A list of the names (including directory) of temporary files created
# for these tests.
filelist = []           # list of names (with directory) of temporary files


def save_file_name(filename):

    global filelist, temp_dir

    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()

    if len(filename) > len(os.path.basename(filename)):
        # `filename` already includes the directory name
        fullname = filename
    else:
        # `filename` does not include the directory yet
        fullname = os.path.join(temp_dir, filename)

    # If the current file is already in `filelist`, don't include it again.
    if fullname not in filelist:
        filelist.append(fullname)

    return fullname


def save_file(fd, filename):

    fullname = save_file_name(filename)

    fd.writeto(fullname, overwrite=True)

    return fullname


def clean_up():

    global filelist, temp_dir

    if temp_dir is None:
        return

    while filelist:
        fullname = filelist.pop()
        os.remove(fullname)
    os.rmdir(temp_dir)
    temp_dir = None


def mk_sci_file(filename, n_imsets=None):

    fd = fits.HDUList(fits.PrimaryHDU())

    dqflags = [4, 8, 512, 0]
    if n_imsets is None:
        n_imsets = 1
    for extver in range(1, n_imsets + 1):
        # SCI HDU
        if extver == 3:
            hdu = fits.ImageHDU(name="SCI", ver=extver)         # no data
        else:
            sci = np.arange(datasize, dtype=np.float32).reshape(shape)
            sci += extver
            hdu = fits.ImageHDU(data=sci.copy(), name="SCI", ver=extver)
        fd.append(hdu)
        # DQ HDU
        if extver == 2:
            hdu = fits.ImageHDU(name="DQ", ver=extver)          # no data
        else:
            i = extver % 4
            dq = np.zeros(shape, dtype=np.int16)
            dq[i, i] = dqflags[i]
            hdu = fits.ImageHDU(data=dq.copy(), name="DQ", ver=extver)
        fd.append(hdu)
        # ERR HDU
        if extver == 1:
            hdu = fits.ImageHDU(name="ERR", ver=extver)         # no data
        else:
            err = np.ones(shape, dtype=np.float32)
            hdu = fits.ImageHDU(data=err.copy(), name="ERR", ver=extver)
        fd.append(hdu)

    fullname = save_file(fd, filename)
    fd.close()

    return fullname


def mk_fringe_flat(filename, hdunum=1):

    if hdunum not in [0, 1]:
        print('Warning (mk_fringe_flat):  hdunum must be 0 or 1')

    if hdunum == 0:
        data = np.ones(shape, dtype=np.float32) * 2.
        fd = fits.HDUList(fits.PrimaryHDU(data=data))
    else:
        extver = 1
        fd = fits.HDUList(fits.PrimaryHDU())
        data = np.ones(shape, dtype=np.float32) * 2.
        data[1, 4] = -1.
        data[1, 5] = 0.
        hdu = fits.ImageHDU(data=data.copy(), name="SCI", ver=extver)
        fd.append(hdu)
        dq = np.zeros(shape, dtype=np.int16)
        dq[2, 2:4] = 15
        hdu = fits.ImageHDU(data=dq.copy(), name="DQ", ver=extver)
        fd.append(hdu)

    fullname = save_file(fd, filename)
    fd.close()

    return fullname


def test_defringe():

    n_imsets = 5

    science_file = mk_sci_file('ttt_crj.fits', n_imsets=n_imsets)

    fringe_flat = mk_fringe_flat('fringe_flat.fits', hdunum=1)

    drj_filename = defringe(science_file, fringe_flat,
                            overwrite=True, verbose=True)
    drj_filename = save_file_name(drj_filename)

    in_sci = fits.open(science_file)
    in_flat = fits.open(fringe_flat)
    out_sci = fits.open(drj_filename)

    fringe_data = in_flat[('sci', 1)].data
    fringe_dq = in_flat[('dq', 1)].data

    results = []

    # Test 1
    check = len(out_sci) == len(in_sci)
    if not check:
        print('FAILED:  the input and output science files do not have the same number of HDUs')
    results.append(check)

    # Test 2
    # Check that the science data array was divided by the fringe flat,
    # except at pixels where the fringe flat was less than or equal to zero,
    # or the fringe flat was flagged as bad in its DQ array.
    f_flat = fringe_data.copy()         # because we're going to modify it
    sci5_data = in_sci[('sci', 5)].data
    ff_is_bad = (fringe_data <= 0.)
    f_flat[ff_is_bad] = 1.
    dq_flagged = (fringe_dq > 0)
    f_flat[dq_flagged] = 1.
    compare = sci5_data / f_flat
    check = np.allclose(out_sci[('sci', 5)].data, compare, rtol=1.e-6)
    if not check:
        print('FAILED:  science array not flat fielded correctly')
    results.append(check)

    # Test 3
    f_flat = fringe_data.copy()
    err4_data = in_sci[('err', 4)].data
    ff_is_bad = (fringe_data <= 0.)
    f_flat[ff_is_bad] = 1.
    dq_flagged = (fringe_dq > 0)
    f_flat[dq_flagged] = 1.
    compare = err4_data / f_flat
    check = np.allclose(out_sci[('err', 4)].data, compare, rtol=1.e-6)
    if not check:
        print('FAILED:  error array not flat fielded correctly')
    results.append(check)

    # Test 4
    # The DQ extension in imset 2 in the input science file did not have a
    # data block, but one should have been created and added to this imset
    # in the output science file.
    dq2_data = out_sci[('dq', 2)].data
    if dq2_data is None:
        print("FAILED:  out_sci[('dq', 2)].data is None")
    results.append(dq2_data is not None)

    # Test 5
    dq4_data = out_sci[('dq', 4)].data
    # Pixels [1, 4], [1, 5], [2, 2], and [2, 3] should all be flagged with 512.
    # The first two were negative or zero in the fringe flat data, and the
    # latter two were flagged with 15 (4 and 8 are both "serious") in the
    # fringe flat data quality array.
    # The input science DQ array had value 4 at pixel [0, 0] for imset 4, and
    # that value should be propagated to the output DQ array.
    compare = np.zeros(shape, dtype=np.int16)
    compare[1, 4] = 512
    compare[1, 5] = 512
    compare[2, 2] = 512 | 15
    compare[2, 3] = 512 | 15
    compare[0, 0] = 4
    diff = dq4_data - compare
    if len(diff.nonzero()[0]) > 0:
        print("FAILED:  out_sci[('dq', 4)].data is not what was expected")
    results.append(len(diff.nonzero()[0]) == 0)

    in_sci.close()
    in_flat.close()
    in_sci.close()

    clean_up()                          # remove temporary files and directory
