#! /usr/bin/env python

import os
import sys
import getopt
import glob
import subprocess

import numpy.random as rn               # used by mkRandomName
from astropy.io import fits

from stsci.tools import parseinput, teal

"""
Perform wavelength calibration of STIS data.

Examples
--------

In Python without TEAL:

>>> import stistools
>>> stistools.wavecal.wavecal("o66p01020_flt.fits", "o66p01020_wav.fits",
...             verbose=True, trailer="o66p01020.trl")

In Python with TEAL:

>>> from stistools import wavecal
>>> from stsci.tools import teal
>>> teal.teal("wavecal")

In Pyraf:

>>> import stistools
>>> teal wavecal

From command line::

% ./wavecal.py -v -s o66p01020_flt.fits o66p01020_wav.fits
% ./wavecal.py -v -s o66p01020_flt.fits o66p01020_w2d_tmp.fits
% ./wavecal.py -r

"""

__taskname__ = "wavecal"
__version__ = "3.4"
__vdate__ = "13-November-2013"
__author__ = "Phil Hodge, STScI, November 2013."

# MJD after which the external shutter was closed for CCD HITM wavecals.
SH_CLOSED = 51126.0


def main(args):

    if len(args) < 2:
        prtOptions()
        print("Specify at least a calibrated science file and its wavecal.")
        sys.exit()

    try:
        (options, pargs) = getopt.getopt(args, "srtv:",
                                         ["version"])
    except Exception as error:
        prtOptions()
        sys.exit()

    input = ""
    inwave = ""
    savetmp = False
    verbose = False
    timestamps = False

    rn.seed()           # used by mkRandomName

    for i in range(len(options)):
        if options[i][0] == "--version":
            status = subprocess.call(["cs4.e", "--version"])
            return 0
        if options[i][0] == "-r":
            status = subprocess.call(["cs4.e", "-r"])
            return 0
        elif options[i][0] == "-v":
            verbose = True
        elif options[i][0] == "-t":
            timestamps = True
        elif options[i][0] == "-s":
            savetmp = True

    nargs = len(pargs)
    if nargs < 1 or nargs > 2:
        prtOptions()
        sys.exit()
    input = pargs[0]
    if nargs == 2:
        outroot = pargs[1]

    status = wavecal(input, wavecal=inwave, debugfile="",
                     savetmp=savetmp,
                     option="linear", angle=None,
                     verbose=verbose, timestamps=timestamps)

    sys.exit(status)


def prtOptions():
    """Print a list of command-line options and arguments."""

    print("The command-line options are:")
    print("  --version (print the version number and exit)")
    print("  -r (print the full version string and exit)")
    print("  -v (verbose)")
    print("  -t (print timestamps)")
    print("  -s (save temporary files)")
    print("")
    print("Following the options, list the input flt file names and")
    print("  the associated raw (or calibrated) wavecal file names.")


def wavecal(input, wavecal, debugfile="", savetmp=False,
            option="linear", angle=None,
            verbose=False, timestamps=False,
            trailer="", print_version=False, print_revision=False):
    """Perform wavecal processing for STIS data.

    Parameters
    ----------
    input: str
        Names of the flt or crj file or files for the science exposure.
        The SHIFTA1 and SHIFTA2 keywords will be updated in these files,
        based on the results of processing the wavecal file(s).

    wavecal: str
        Names of the associated wavecal file or files (either raw or
        calibrated).  If this is a raw file, it will first be calibrated
        using cs1.e (basic2d), then with cs7.e (x2d) for first-order
        grating data.  These calibrated files are regarded as temporary,
        and (unless savetmp) they will be deleted when processing has
        been completed.

    debugfile: str
        If specified, debugging information will be written to a file with
        this name.  For echelle data this will be a FITS file, but for
        first-order data it will be a text file (and possibly a FITS file
        as well).

    savetmp: bool
        If wavecal is a raw wavecal file, some calibration will be
        performed, depending on mode.  If savetmp is False (the default),
        the calibrated wavecal files will be deleted after wavecal
        processing is complete.

    option: str
        If the wavecal file contains more than one image set, the shifts
        will be interpolated between wavecal exposures that bracket the
        science exposure.  This argument gives the interpolation option,
        either "linear" (the default) or "nearest".  If the science
        exposure was before the first or after the last exposure in the
        wavecal file, the shifts will be copied from the first or last
        exposure respectively.

    angle: float or None
        This argument is only relevant for echelle data for which the
        wavecal was taken with a long slit (e.g. 6X0.2).  The angles have
        not been measured accurately; they vary from one grating to
        another, and they even vary depending on location on the detector.
        This argument specifies the slit angle, in degrees measured
        clockwise from the Y axis.  Here are some approximate values:

            - E230M:  0.9 to 1.2
            - E230H:  4.9 to 6.9
            - E140H: -3.8 to -5.8

    verbose: bool
        If True, calstis will print more info.

    timestamps: bool
        If True, calstis will print the date and time at various points
        during processing.

    trailer: str
        If specified, the standard output and standard error will be
        written to this file instead of to the terminal.  Note, however,
        that if print_version or print_revision is specified, the value
        will be printed to the terminal, and any name given for the
        trailer will be ignored.

    print_version: bool
        If True, calstis will print the version number (a string) and
        then return 0, without checking any other argument.

    print_revision: bool
        If True, calstis will print the full version string and then
        return 0.

    Returns
    -------
    status: int
        0 is OK.
        1 is returned if cs4.e (the calstis host executable) returned a
        non-zero status.  If verbose is True, the value returned by cs4.e
        will be printed.
        2 is returned if the specified input file or files were not found,
        of if the numbers of input and wavecal files (or of debugfiles) are
        not the same.
    """

    if print_version:
        status = subprocess.call(["cs4.e", "--version"])
        return 0
    if print_revision:
        status = subprocess.call(["cs4.e", "-r"])
        return 0

    # infiles may include one or more file names, separated by blanks
    # or commas (or both), and any name may include wildcards.
    infiles = []
    input1 = input.split()
    for in1 in input1:
        input2 = in1.split(",")
        for in2 in input2:
            files = glob.glob(in2)
            infiles.extend(files)
    if input1 and not infiles:
        print("No file name matched the string '{}'".format(input))
        return 2

    wavecal_files = []
    wavecal1 = wavecal.split()
    for wav1 in wavecal1:
        if wav1:
            wavecal2 = wav1.split(",")
            for wav2 in wavecal2:
                if wav2:
                    wavecal_files.append(wav2)

    dbgfiles = []
    if debugfile:
        dbgfiles1 = debugfile.split()
        for dbg1 in dbgfiles1:
            dbgfiles2 = dbg1.split(",")
            for dbg2 in dbgfiles2:
                dbgfiles.append(dbg2)

    same_length = True          # optimistic initial value
    n_infiles = len(infiles)
    if wavecal_files and len(wavecal_files) != n_infiles:
        same_length = False
        print("You specified {} input files but {} wavecal files.".format
              (n_infiles, len(wavecal_files)))
        print("The number of input and wavecal files must be the same.")
    if dbgfiles and len(dbgfiles) != n_infiles:
        same_length = False
        print("The number of input and debugfile files must be the same.")
    if not same_length:
        return 2

    if trailer:
        if verbose and os.access(trailer, os.F_OK):
            print("Appending to trailer file {}".format(trailer))
        f_trailer = open(trailer, "a")
        fd_trailer = f_trailer.fileno()
    else:
        f_trailer = None
        fd_trailer = None

    for (i, infile) in enumerate(infiles):

        tempfnames = []
        # Run cs1.e on the wavecal file.
        (flag, fwv_file) = runBasic2d(wavecal_files[i], tempfnames,
                                      verbose, timestamps, fd_trailer)

        # Subtract a fraction of the infile (flt or crj) from the wavecal,
        # if the exposure was taken long enough ago.
        (flag, cwv_file) = runCs11(fwv_file, infile, tempfnames,
                                   verbose, timestamps, fd_trailer)

        # Run cs7.e on the wavecal flt file (except for echelle or prism).
        (flag, w2d_file) = runX2d(cwv_file, angle, tempfnames,
                                  verbose, timestamps, fd_trailer)

        # Now run cs4.e on the w2d_file to find the shifts.
        if dbgfiles:
            dbg = dbgfiles[i]
        else:
            dbg = None
        runWavecal(w2d_file, dbg, angle, verbose, timestamps, fd_trailer)

        # Run cs12.e to copy the shifts to infile.
        runCs12(w2d_file, infile, option, verbose, timestamps, fd_trailer)

        if not savetmp:
            for tmp_file in tempfnames:
                if verbose:
                    print(" ... deleting temporary file {}".format(tmp_file))
                try:
                    os.remove(tmp_file)
                except OSError:
                    print("Warning:  couldn't delete temporary file {}.".format
                          (tmp_file))

    if f_trailer is not None:
        f_trailer.close()

    return 0


def mkRandomNameW(prefix="wavecal_", suffix="_tmp.fits", n=100000000):
    MAX_TRIES = 100
    done = False
    k = 0
    while not done:
        i = rn.randint(0, n, 1)[0]
        filename = "%s%d%s" % (prefix, i, suffix)
        k += 1
        if not os.access(filename, os.F_OK):
            done = True
        if k > MAX_TRIES:
            break

    if done:
        return filename
    else:
        return None


def runBasic2d(wavecal, tempfnames, verbose, timestamps, fd_trailer):

    flag = False                # initial value

    # First check whether the wavecal file is already calibrated.
    fd = fits.open(wavecal)
    dqicorr = fd[0].header.get("dqicorr", default="missing")
    blevcorr = fd[0].header.get("blevcorr", default="missing")
    darkcorr = fd[0].header.get("darkcorr", default="missing")
    flatcorr = fd[0].header.get("flatcorr", default="missing")
    detector = fd[0].header.get("detector", default="missing")
    fd.close()
    if dqicorr == "COMPLETE" or blevcorr == "COMPLETE" or \
       darkcorr == "COMPLETE" or flatcorr == "COMPLETE":

        # wavecal is already calibrated.
        fwv_file = wavecal
        flag = False

    else:

        # Create pseudo-random fwv_file name.
        prefix = "wavecal_"
        suffix = "_fwv_tmp.fits"
        fwv_file = mkRandomNameW(prefix, suffix)
        if fwv_file is None:
            raise RuntimeError("Couldn't create temp file name"
                               " %s<digits>%s" % (prefix, suffix))

        arglist = ["cs1.e"]
        if verbose:
            arglist.append("-v")
        if timestamps:
            arglist.append("-t")

        arglist.append(wavecal)
        arglist.append(fwv_file)

        arglist.append("-dqi")
        if detector == "CCD":
            arglist.append("-blev")
            arglist.append("-bias")
        else:
            arglist.append("-lors")
        arglist.append("-dark")
        arglist.append("-flat")

        if verbose:
            print("Running cs1.e on {}".format(wavecal))
            print("  {}".format(arglist))
        status = subprocess.call(arglist, stdout=fd_trailer,
                                 stderr=subprocess.STDOUT)
        if status:
            raise RuntimeError("status = {} from cs1.e".format(status))
        tempfnames.append(fwv_file)
        flag = True

    return flag, fwv_file


def runCs11(fwv_file, infile, tempfnames, verbose, timestamps, fd_trailer):
    """Subtract a fraction of the science image from the wavecal image."""

    # Check whether we need to run cs11.e.
    fd = fits.open(fwv_file)
    sclamp = fd[0].header.get("sclamp", default="missing")
    detector = fd[0].header.get("detector", default="missing")
    texpstrt = fd[0].header.get("texpstrt", default="missing")
    fd.close()
    if detector == "CCD" and sclamp.startswith("HITM") and \
       texpstrt <= SH_CLOSED:
        # Create pseudo-random cwv_file name.
        prefix = "wavecal_"
        suffix = "_cwv_tmp.fits"
        cwv_file = mkRandomNameW(prefix, suffix)
        if cwv_file is None:
            raise RuntimeError("Couldn't create temp file name"
                               " %s<digits>%s" % (prefix, suffix))
        arglist = ["cs11.e"]
        if verbose:
            arglist.append("-v")
        if timestamps:
            arglist.append("-t")

        arglist.append(fwv_file)
        arglist.append(infile)
        arglist.append(cwv_file)

        if verbose:
            print("Running cs11.e")
            print("  {}".format(arglist))
        status = subprocess.call(arglist, stdout=fd_trailer,
                                 stderr=subprocess.STDOUT)
        if status:
            raise RuntimeError("status = {} from cs11.e".format(status))
        tempfnames.append(cwv_file)
        flag = True
    else:
        cwv_file = fwv_file
        flag = False

    return flag, cwv_file


def runX2d(cwv_file, angle, tempfnames,
           verbose, timestamps, fd_trailer):

    flag = False                # initial value

    fd = fits.open(cwv_file)
    opt_elem = fd[0].header.get("opt_elem", default="missing")
    x2dcorr = fd[0].header.get("x2dcorr", default="missing")
    fd.close()

    # Skip 2-D rectification for echelle or prism data.
    skip_it = opt_elem.startswith("E") or opt_elem == "PRISM"

    if skip_it:
        w2d_file = cwv_file
        flag = False
    elif x2dcorr == "COMPLETE":
        # The wavecal is already fully calibrated.
        w2d_file = cwv_file
        flag = False
    else:
        # Create pseudo-random w2d_file name.
        prefix = "wavecal_"
        suffix = "_w2d_tmp.fits"
        w2d_file = mkRandomNameW(prefix, suffix)
        if w2d_file is None:
            raise RuntimeError("Couldn't create temp file name"
                               " %s<digits>%s" % (prefix, suffix))

        arglist = ["cs7.e"]
        if verbose:
            arglist.append("-v")
        if timestamps:
            arglist.append("-t")
        # Specify -x2d and no other calibration switch so that fluxcorr
        # will not be done.
        arglist.append("-x2d")
        arglist.append(cwv_file)
        arglist.append(w2d_file)
        if angle is not None:
            arglist.append("-angle")
            arglist.append("%.20g" % angle)

        if verbose:
            print("Running x2d on {}".format(cwv_file))
            print("  {}".format(arglist))
        status = subprocess.call(arglist, stdout=fd_trailer,
                                 stderr=subprocess.STDOUT)
        if status:
            raise RuntimeError("status = %d from cs7.e" % status)
        tempfnames.append(w2d_file)
        flag = True

    return flag, w2d_file


def runWavecal(w2d_file, dbg, angle, verbose, timestamps, fd_trailer):

    arglist = ["cs4.e"]
    if verbose:
        arglist.append("-v")
    if timestamps:
        arglist.append("-t")
    arglist.append(w2d_file)
    if angle is not None:
        arglist.append("-angle")
        arglist.append("%.20g" % angle)
    if dbg:                             # text file for debug output
        arglist.append("-d")
        arglist.append(dbg)

    if verbose:
        print("Running cs4.e on {}".format(w2d_file))
        print("  {}".format(arglist))
    status = subprocess.call(arglist, stdout=fd_trailer,
                             stderr=subprocess.STDOUT)
    if status:
        raise RuntimeError("status = %d from cs4.e" % status)


def runCs12(w2d_file, infile, option, verbose, timestamps, fd_trailer):

    arglist = ["cs12.e"]
    if verbose:
        arglist.append("-v")
    if timestamps:
        arglist.append("-t")

    arglist.append(w2d_file)
    arglist.append(infile)
    arglist.append(option)

    if verbose:
        print("Running cs12.e")
        print("  {}".format(arglist))
    status = subprocess.call(arglist, stdout=fd_trailer,
                             stderr=subprocess.STDOUT)
    if status:
        raise RuntimeError("status = %d from cs12.e" % status)

#-------------------------#
# Interfaces used by TEAL #
#-------------------------#


def getHelpAsString(fulldoc=True):
    """Return documentation on the wavecal function."""
    return wavecal.__doc__


def run(configobj=None):
    """TEAL interface for the wavecal function."""
    wavecal(input=configobj["input"],
            wavecal=configobj["wavecal"],
            debugfile=configobj["debugfile"],
            savetmp=configobj["savetmp"],
            option=configobj["option"],
            angle=configobj["angle"],
            verbose=configobj["verbose"],
            timestamps=configobj["timestamps"],
            trailer=configobj["trailer"],
            print_version=configobj["print_version"],
            print_revision=configobj["print_revision"])

if __name__ == "__main__":

    main(sys.argv[1:])
