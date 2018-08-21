#!/usr/bin/env python
"""
A Python module for aligning the spectra in different flat-fielded
images of an IMSET.  These files can then be combined with
along-the-slit dithering to reject hot pixels and cosmic rays.  The
POSTARG2 keyword is used to determine the number of rows to be
shifted.
"""

from astropy.io import fits


__version__ = '1.7 (2010-Apr-27)'


def shiftimage(infile, outfile, shift=0):

    """
    Shift each image extension of an input file by N rows and write the
    new image extension to the output file.
    """

    fin = fits.open(infile)           # flat-fielded file

    fout = fits.HDUList()              # shifted flat-field file
    phdr = fin[0].header
    phdr.add_history('SSHIFT complete ...')
    phdr.add_history('  all extensions were shifted by %d rows' % shift)
    fout.append(fits.PrimaryHDU(header=phdr))

    for exten in fin[1:]:
        image = exten.data.copy()
        image[:, :] = 0
        if shift > 0:
            image[shift:] = exten.data[:-shift]
        elif shift < 0:
            image[:shift] = exten.data[-shift:]
        else:
            image[:] = exten.data[:]
        fout.append(fits.ImageHDU(header=exten.header, data=image))

    fout.writeto(outfile)


def sshift(input, output=None, shifts=None, platescale=None,
           tolerance=None):

    """ Align spectra from different images of an imset.

Parameters
----------
input : list
    A list of input filenames.  These must be STIS flat-
    fielded (_flt) image FITS files.  This argument will accept a
    single filename or a list of filenames.
shifts : list, optional
    A list of integers indicating the number of rows to shift
    each image of each file in the cross-dispersion (Y-) direction.
platescale : float, optional
    The size of a pixel in arcseconds.  Used to convert
    the value of the POSTARG2 keyword to pixels.
tolerance : float, optional
    The allowed difference between calculated shifts and
    integer pixel shifts (fraction of pixel).

Returns
-------
output : list, optional
    A list of output filenames. The number of output
    filenames must match the number of input filenames.  If no output
    is given, then the _flt substring of the input file is replace by
    the _sfl substring to create an output file.  This option will
    accept a single filename or a list of filenames.

Notes
------
Author:
  - Paul Barrett (STScI)

    """

    # History:
    # 2003/09/22  PEB - version 1.0
    # 2003/11/05  PEB - version 1.1
    # 2003/11/10  PEB - version 1.2 - add postarg1/2, wavecal checks
    # 2003/11/17  PEB - version 1.3 - add history cards
    # 2004/09/13  PEB - version 1.4 - set PLATESC, default:0.0507
    #                                 add tolerance keyword, default: 0.1
    # 2004/09/23  PEB - version 1.5 - check for mixed dataset
    #                               - fixed integer shift bug
    #                                 removed wavecorr step
    #                                 check for binned data and non-integral
    #                                 shifts.
    # 2004/09/24  PEB - version 1.6 - add keyword consistency checks
    #  Setup input and output filename lists, so iteration can be done
    #  over a list of zipped filenames.

    if not isinstance(input, list):
        input = [input]
    elif not input:
        raise ValueError(
            'No input files found.  Possibly using wrong directory.')

    if output is None:
        output = len(input)*[None]
    elif not isinstance(output, list):
        output = [output]

    if shifts is None:
        pass
    elif not isinstance(shifts, list):
        shifts = [shifts]

    if tolerance is None:
        tolerance = 0.1

    if len(input) != len(output):
        raise ValueError(
            'number of output files is not equal to number input files')

    if shifts is not None:
        for shift in shifts:
            if not isinstance(shift, int):
                raise ValueError('shift value must be an integer')

    xposs, yposs, xpos0, ypos0 = [], [], None, None
    proposid, obset_id, targname = None, None, None
    propaper, opt_elem, cenwave = None, None, None
    binaxis1, binaxis2 = None, None
    for infile in input:

        #  Read the POSTARG2 keyword (in the primary header) and the
        #  CRPIX2 keyword (in the 1st extension) to determine the
        #  relative shift of each input file.  Choose a reference
        #  position that is closest to Y-pixel 512.

        infil = fits.open(infile)
        phdr = infil[0].header

        if platescale is None:
            # platescale = phdr['PLATESC']
            platescale = 0.05077
        else:
            platescale = float(platescale)

        if phdr['FLATCORR'].upper() != 'COMPLETE':
            raise ValueError(
                'Input file has not been flat-fielded corrected.')

        #  Check that TARGNAME is the same.
        if targname is None:
            targname = phdr['TARGNAME']
        elif targname != phdr['TARGNAME']:
            raise ValueError('Not all exposures are for the same target.')

        #  Check that all PROPOSID and OBSET_ID values are the same.
        if proposid is None:
            proposid = phdr['PROPOSID']
            obset_id = phdr['OBSET_ID']
        elif proposid != phdr['PROPOSID'] or obset_id != phdr['OBSET_ID']:
            raise ValueError(' Not all exposures are from the same visit;'
                             ' placement of the spectrum on the detector will'
                             ' differ.')

        # Check that PROPAPER, OPT_ELEM, CENWAVE are the same.
        if propaper is None:
            propaper = phdr['PROPAPER']
            opt_elem = phdr['OPT_ELEM']
            cenwave = phdr['CENWAVE']
        elif propaper != phdr['PROPAPER'] or opt_elem != phdr['OPT_ELEM'] or \
             cenwave != phdr['CENWAVE']:
            raise ValueError('Different observing configurations have been used.')

        #  Check that BINAXIS1 and BINAXIS2 are the same.
        if binaxis1 is None:
            binaxis1 = phdr['BINAXIS1']
            binaxis2 = phdr['BINAXIS2']
        elif binaxis1 != phdr['BINAXIS1'] or binaxis2 != phdr['BINAXIS2']:
            raise ValueError('Different binnings have been used.')

        #  Check that all POSTARG1 values are the same (within reason).
        xpos = phdr['POSTARG1']
        if xpos0 is None:
            xpos0 = xpos
        elif abs(xpos - xpos0) > 0.05:
            raise ValueError('POSTARG1 values of input files are not equal.')

        #  Get the POSTARG2 values and the one that is nearest to row 512.
        ypos = phdr['POSTARG2'] / platescale
        ypix = infil[1].header['CRPIX2'] - 512
        if ypos0 is None or abs(ypix + ypos) < abs(ypix + ypos0):
            ypos0 = ypos

        yposs.append(ypos)

    #  Check for non-integral POSTARG2 values and calculate array of
    #  pixel shifts.
    if shifts is None:
        shifts = []
        for ypos in yposs:
            dypos = ypos - ypos0
            if abs(abs(dypos) - int(abs(dypos)+0.5)) > tolerance:
                raise ValueError("POSTARG2 shift not within the specified "
                                 "tolerance {} pix of integer-pixel shift".format(tolerance))
            # 'POSTARG2 shift greater than specified tolerance: %d' % tolerance
            # 'non-integral POSTARG2 value or incorrect plate scale.'
            if dypos < 0.:
                ishift = -int(dypos-0.5)
            else:
                ishift = -int(dypos+0.5)
            if ishift % binaxis2:
                raise ValueError('Non-integral pixel shift for binned data')
            shifts.append(ishift//binaxis2)

    #  Process each file using corresponding pixel shift.
    print('input-file        pixel-shift')
    for infile, outfile, npixel in zip(input, output, shifts):
        fin = fits.open(infile)

        #  Use default output file name.
        if outfile is None:
            import re
            outfile = re.sub('flt\.', 'sfl.', infile, count=1)

        if binaxis2 == 1:
            print('{:>18}: {:3}'.format(infile, npixel))
        else:
            print('{:>18}: {:3}  binned'.format(infile, npixel))
        shiftimage(infile, outfile, shift=npixel)
        fin.close()


if __name__ == '__main__':
    import sys
    import getopt

    output, shifts, scale, toler = None, None, None, None

    short_opts = 'o:s:p:t:h'
    long_opts = ['output=', 'shifts=', 'platescale=', 'tolerance=', 'help']

    opts, args = getopt.getopt(sys.argv[1:], short_opts, long_opts)

    for opt in opts:
        if opt[0] == '-o' or opt[0] == '--output':
            output = opt[1].split(',')
        elif opt[0] == '-s' or opt[0] == '--shifts':
            shifts = [int(v) for v in opt[1].split(',')]
        elif opt[0] == '-p' or opt[0] == '--platescale':
            scale = eval(opt[1])
        elif opt[0] == '-t' or opt[0] == '--tolerance':
            toler = eval(opt[1])
        elif opt[0] == '-h' or opt[0] == '--help':
            print(sshift.__doc__)
            sys.exit()

    if len(args) > 0:
        sshift(args, output=output, shifts=shifts, platescale=scale,
               tolerance=toler)
    else:
        print("""Usage: sshift [-o|--output 'files'] [-s|--shifts 'shifts']
        [-p|--platescale scale] [-t|--tolerance tol] [-h|--help] input-files""")
