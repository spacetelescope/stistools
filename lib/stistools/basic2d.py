#! /usr/bin/env python

"""
Calibrate STIS data.

Examples
--------

In Python without TEAL:

>>> import stistools
>>> stistools.basic2d.basic2d("o66p01020_raw.fits", verbose=True,
...                           trailer="o66p01020.trl")

In Python with TEAL:

>>> from stistools import basic2d
>>> from stsci.tools import teal
>>> teal.teal("basic2d")

In Pyraf:

>>> import stistools
>>> teal basic2d

From command line::

% ./basic2d.py -v -s o66p01020_raw.fits o66p01020_flt.fits
% ./basic2d.py -r
"""

from __future__ import division         # confidence unknown
import os
import sys
import getopt
import glob
import subprocess

from stsci.tools import parseinput,teal

__taskname__ = "basic2d"
__version__ = "3.0"
__vdate__ = "14-January-2013"
__author__ = "Phil Hodge, STScI, January 2013."

def main(args):

    if len(args) < 1:
        prtOptions()
        print("At least a raw file name must be specified.")
        sys.exit()

    try:
        (options, pargs) = getopt.getopt(args, "rtv:",
                                         ["version"])
    except Exception, error:
        prtOptions()
        sys.exit()

    output = ""
    outblev = ""
    verbose = False
    timestamps = False

    for i in range(len(options)):
        if options[i][0] == "--version":
            status = subprocess.call(["cs1.e", "--version"])
            return 0
        if options[i][0] == "-r":
            status = subprocess.call(["cs1.e", "-r"])
            return 0
        elif options[i][0] == "-v":
            verbose = True
        elif options[i][0] == "-t":
            timestamps = True

    nargs = len(pargs)
    if nargs < 1 or nargs > 2:
        prtOptions()
        sys.exit()
    input = pargs[0]
    if nargs == 2:
        output = pargs[1]

    status = basic2d(input, output=output, outblev=outblev,
                     verbose=verbose, timestamps=timestamps)

    sys.exit(status)

def prtOptions():
    """Print a list of command-line options and arguments."""

    print("The command-line options are:")
    print("  --version (print the version number and exit)")
    print("  -r (print the full version string and exit)")
    print("  -v (verbose)")
    print("  -t (print timestamps)")
    print("")
    print("Following the options, list one or more input raw file names,")
    print("  enclosed in quotes if more than one file name is specified")
    print("  and/or if wildcards are used.")
    print("One or more output file names may be specified (the same number")
    print("  as the input file names).")

def basic2d(input, output="", outblev="",
           dqicorr=True, blevcorr=True, doppcorr=True,
           lorscorr=True, glincorr=True, lflgcorr=True,
           biascorr=True, darkcorr=True, flatcorr=True,
           photcorr=True, statflag=True,
           darkscale="",
           verbose=False, timestamps=False,
           trailer="", print_version=False, print_revision=False):
    """Perform basic 2-D calibration of STIS raw data.

    Some calibration steps are relevant only for CCD or only for MAMA, and
    since an output file of calstis or basic2d may be used as the input,
    some steps may have already been done.  Most calibration steps will not
    be done if they are not relevant or if they have already been done,
    regardless of the value of the calibration switch (e.g. flatcorr).

    Parameters
    ----------
    input: str
        Name of the input raw file.

    output: str
        Name of the output file, or "" (the default).  If no name was
        specified, the output name will be constructed from the input name.

    outblev: str
        Name of the output text file for blev info, or "" (the default).

    dqicorr: bool
        If True, update the DQ array.

    blevcorr: bool
        If True, subtract a bias level based on the overscan values.
        (CCD only.)

    doppcorr: bool
        If True, convolve reference files (bpixtab, darkfile, flatfile)
        as needed with the Doppler shift offset throughout the exposure,
        if Doppler correction was done on-board.  (MAMA only, because
        for the CCD Doppler correction is not done on-board.)

    lorscorr: bool
        If True, bin high-res data to lo-res.  (MAMA only.)

    glincorr: bool
        If True, correct for global non-linearity.  (MAMA only.)

    lflgcorr: bool
        If True, flag local non-linearity.  (MAMA only.)

    biascorr: bool
        If True, subtract the bias image.  (CCD only.)

    darkcorr: bool
        If True, subtract the dark image, scaled by the exposure time
        and possibly also a temperature-dependent factor.

    flatcorr: bool
        If True, divide by the flat field image.

    photcorr: bool
        If True, determine the photometric parameters and populate keywords
        PHOTFLAM, PHOTZPT, PHOTPLAM and PHOTBW.  (Imaging data only.)

    statflag: bool
        If True, compute statistics for image arrays and update keywords.

    darkscale: str
        This may be used to override the time and/or temperature dependent
        scale factor that would normally be applied to the dark image
        before subtracting from the raw data.  It's a string rather than
        a float in order to accept a different scale factor for each
        image set in the input data.  calstis reads the value or values
        (separated by blanks) from the string, and if the value is greater
        than zero, it will be used instead of the value determined from
        the temperature and time.  (CCD or NUV-MAMA only.)

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
        then return 0.

    print_revision: bool
        If True, calstis will print the full version string and then
        return 0.

    Returns
    -------
    status: int
        0 is OK.
        1 is returned if cs1.e (the calstis host executable) returned a
        non-zero status.  If verbose is True, the value returned by cs1.e
        will be printed.
        2 is returned if the specified input file or files were not found,
        or if there is a mismatch between the number of input, output,
        and/or outblev files specified.
    """

    if print_version:
        status = subprocess.call(["cs1.e", "--version"])
        return 0
    if print_revision:
        status = subprocess.call(["cs1.e", "-r"])
        return 0

    cumulative_status = 0

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
        print("No file name matched the string '%s'" % input)
        return 2

    if output:
        outfiles = []
        output1 = output.split()
        for out1 in output1:
            output2 = out1.split(",")
            for out2 in output2:
                outfiles.append(out2)
    else:
        outfiles = None

    if outblev:
        outblev_txt = []
        outblev1 = outblev.split()
        for out1 in outblev1:
            outblev2 = out1.split(",")
            for out2 in outblev2:
                outblev_txt.append(out2)
    else:
        outblev_txt = None

    same_length = True          # optimistic initial value
    n_infiles = len(infiles)
    if outfiles and len(outfiles) != n_infiles:
        same_length = False
        print("You specified %d input files but %d output files." %
              (n_infiles, len(outfiles)))
        print("The number of input and output files must be the same.")
    if outblev_txt and len(outblev_txt) != n_infiles:
        same_length = False
        print("The number of input and outblev files must be the same.")
    if not same_length:
        return 2

    if trailer is None:
        f_trailer = None
    else:
        f_trailer = open(trailer, "w")

    for (i, infile) in enumerate(infiles):

        arglist = ["cs1.e"]
        if verbose:
            arglist.append("-v")
        if timestamps:
            arglist.append("-t")
        if darkscale:
            arglist.append("-dscl")
            arglist.append("%s" % darkscale)

        arglist.append(infile)
        if outfiles is not None:
            arglist.append(outfiles[i])
        if outblev_txt is not None:
            arglist.append(outblev_txt[i])

        if dqicorr:
            arglist.append("-dqi")
        if blevcorr:
            arglist.append("-blev")
        if doppcorr:
            arglist.append("-dopp")
        if lorscorr:
            arglist.append("-lors")
        if glincorr:
            arglist.append("-glin")
        if lflgcorr:
            arglist.append("-lflg")
        if biascorr:
            arglist.append("-bias")
        if darkcorr:
            arglist.append("-dark")
        if flatcorr:
            arglist.append("-flat")
        if photcorr:
            arglist.append("-phot")
        if statflag:
            arglist.append("-stat")

        if verbose:
            print("'%s'" % str(arglist))
            print("Running basic2d on %s" % infile)
        if f_trailer is None:           # no trailer file
            status = subprocess.call(arglist)
        else:
            status = subprocess.call(arglist, stdout=f_trailer.fileno(),
                                     stderr=subprocess.STDOUT)
            if status:
                cumulative_status = 1
                if verbose:
                    print("Warning:  status = %d" % status)

    if f_trailer is not None:
        f_trailer.close()

    return cumulative_status

#-------------------------#
# Interfaces used by TEAL #
#-------------------------#

def getHelpAsString(fulldoc=True):
    """Return documentation on the basic2d function."""
    return basic2d.__doc__

def run(configobj=None):
    """TEAL interface for the basic2d function."""
    basic2d(configobj["input"],
            configobj["output"],
            configobj["outblev"],
            configobj["dqicorr"],
            configobj["blevcorr"],
            configobj["doppcorr"],
            configobj["lorscorr"],
            configobj["glincorr"],
            configobj["lflgcorr"],
            configobj["biascorr"],
            configobj["darkcorr"],
            configobj["flatcorr"],
            configobj["photcorr"],
            configobj["statflag"],
            configobj["darkscale"],
            configobj["verbose"],
            configobj["timestamps"],
            configobj["trailer"],
            configobj["print_version"],
            configobj["print_revision"])

if __name__ == "__main__":

    main(sys.argv[1:])
