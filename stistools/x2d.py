#! /usr/bin/env python

import os
import sys
import getopt
import glob
import subprocess

from stsci.tools import parseinput, teal

__doc__ = """
Rectify 2-D STIS spectral data.

Examples
--------

In Python without TEAL:

>>> import stistools
>>> stistools.x2d.x2d("o66p01020_flt.fits", output="test_x2d.fits",
...                   verbose=True, trailer="o66p01020.trl")

In Python with TEAL:

>>> from stistools import x2d
>>> from stsci.tools import teal
>>> teal.teal("x2d")

From command line::

% ./x2d.py -v o66p01020_flt.fits o66p01020_x2d.fits
% ./x2d.py -r
"""

__taskname__ = "x2d"
__version__ = "3.4"
__vdate__ = "13-November-2013"
__author__ = "Phil Hodge, STScI, November 2013."


def main(args):

    if len(args) < 1:
        prtOptions()
        print("At least an input file name must be specified.")
        sys.exit()

    try:
        (options, pargs) = getopt.getopt(args, "rtv:",
                                         ["version"])
    except Exception as error:
        prtOptions()
        sys.exit()

    output = ""
    blazeshift = None
    verbose = False
    timestamps = False

    for i in range(len(options)):
        if options[i][0] == "--version":
            status = subprocess.call(["cs7.e", "--version"])
            return 0
        if options[i][0] == "-r":
            status = subprocess.call(["cs7.e", "-r"])
            return 0
        if options[i][0] == "-v":
            verbose = True
        if options[i][0] == "-t":
            timestamps = True

    nargs = len(pargs)
    if nargs < 1 or nargs > 2:
        prtOptions()
        sys.exit()
    input = pargs[0]
    if nargs == 2:
        output = pargs[1]

    status = x2d(input, output=output,
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
    print("Following the options, list one or more input file names,")
    print("  enclosed in quotes if more than one file name is specified")
    print("  and/or if wildcards are used.")
    print("One or more output file names may be specified (the same number")
    print("  as the input file names).")


def x2d(input, output="",
        helcorr="perform", fluxcorr="perform", statflag=True,
        center=False, blazeshift=None, err_alg="wgt_var",
        verbose=False, timestamps=False, trailer="",
        print_version=False, print_revision=False):
    """Rectify 2-D STIS spectral data.

    Parameters
    ----------
    input: str
        Name of the input raw file.

    output: str
        Name of the output file, or "" (the default).  If no name was
        specified, the output name will be constructed from the input name.

    helcorr: str
        If "perform", correct for heliocentric Doppler shift.

    fluxcorr: str
        If "perform", convert to absolute flux.

    statflag: bool
        If True, compute statistics for image arrays and update keywords.

    center: bool
        If True, center the target spectrum in the cross-dispersion
        direction.  For G140L and G140M spectra, the target has at
        different times been offset to a location either above or below
        the middle of the detector, to avoid the repeller wire.  This
        argument allows more convenient comparison of data taken at
        widely different times.

    blazeshift: float or None
        Blaze shift (in pixels).  None means blazeshift is not specified.

    err_alg: str
        Algorithm for computing error estimates.  The default is "wgt_var",
        which means that the weight (for bilinear interpolation) is applied
        to the variances of the input pixels.  The alternative is
        "wgt_err", to specify that the weight should be applied to the
        errors of the input pixels.

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
        1 is returned if cs7.e (the calstis host executable) returned a
        non-zero status.  If verbose is True, the value returned by cs7.e
        will be printed.
        2 is returned if the specified input file or files were not found,
        or the numbers of input and output files (if the latter was
        specified) are not the same.
    """

    if print_version:
        status = subprocess.call(["cs7.e", "--version"])
        return 0
    if print_revision:
        status = subprocess.call(["cs7.e", "-r"])
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
        print("No file name matched the string '{}'".format(input))
        return 2

    if output:
        outfiles = []
        output1 = output.split()
        for out1 in output1:
            if out1:
                output2 = out1.split(",")
                for out2 in output2:
                    if out2:
                        outfiles.append(out2)
    else:
        outfiles = None

    n_infiles = len(infiles)
    if outfiles and len(outfiles) != n_infiles:
        print("You specified {} input files but {} output files.".format
              (n_infiles, len(outfiles)))
        print("The number of input and output files must be the same.")
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

        arglist = ["cs7.e"]

        arglist.append(infile)
        if outfiles is not None:
            arglist.append(outfiles[i])

        if verbose:
            arglist.append("-v")
        if timestamps:
            arglist.append("-t")

        switch_was_set = False
        if helcorr == "perform":
            arglist.append("-hel")
            switch_was_set = True
        if fluxcorr == "perform":
            arglist.append("-flux")
            switch_was_set = True
        if not switch_was_set:
            arglist.append("-x2d")

        if err_alg:
            if err_alg == "wgt_err":
                arglist.append("-wgt_err")
            elif err_alg != "wgt_var":
                raise RuntimeError("err_alg must be either 'wgt_err'"
                                   " or 'wgt_var'; you specified '%s'" % err_alg)

        if blazeshift is not None:
            arglist.append("-b")
            arglist.append("%.10g" % blazeshift)

        if verbose:
            print("Running x2d on {}".format(infile))
            print("  {}".format(arglist))
        status = subprocess.call(arglist, stdout=fd_trailer,
                                 stderr=subprocess.STDOUT)
        if status:
            cumulative_status = 1
            if verbose:
                print("Warning:  status = {}".format(status))

    if f_trailer is not None:
        f_trailer.close()

    return cumulative_status

#-------------------------#
# Interfaces used by TEAL #
#-------------------------#


def getHelpAsString(fulldoc=True):
    """Return documentation on the x2d function."""
    return x2d.__doc__


def run(configobj=None):
    """TEAL interface for the x2d function."""
    x2d(input=configobj["input"],
        output=configobj["output"],
        helcorr=configobj["helcorr"],
        fluxcorr=configobj["fluxcorr"],
        statflag=configobj["statflag"],
        center=configobj["center"],
        blazeshift=configobj["blazeshift"],
        err_alg=configobj["err_alg"],
        verbose=configobj["verbose"],
        timestamps=configobj["timestamps"],
        trailer=configobj["trailer"],
        print_version=configobj["print_version"],
        print_revision=configobj["print_revision"])

if __name__ == "__main__":

    main(sys.argv[1:])
