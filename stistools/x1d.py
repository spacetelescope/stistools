#! /usr/bin/env python

import os
import sys
import getopt
import glob
import subprocess

from stsci.tools import parseinput, teal

__doc__ = """
Extract 1-D spectrum.

Examples
--------

In Python without TEAL:

>>> import stistools
>>> stistools.x1d.x1d("o66p01020_flt.fits", output="test_x1d.fits",
...                   verbose=True, trailer="o66p01020.trl")

In Python with TEAL:

>>> from stistools import x1d
>>> from stsci.tools import teal
>>> teal.teal("x1d")

From command line::

% ./x1d.py -v o66p01020_flt.fits o66p01020_x1d.fits
% ./x1d.py -r
"""

__taskname__ = "x1d"
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
    verbose = False
    timestamps = False

    for i in range(len(options)):
        if options[i][0] == "--version":
            status = subprocess.call(["cs6.e", "--version"])
            return 0
        if options[i][0] == "-r":
            status = subprocess.call(["cs6.e", "-r"])
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

    status = x1d(input, output=output,
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
    print("Following the options, list one or more flt or crj file names,")
    print("  enclosed in quotes if more than one file name is specified")
    print("  and/or if wildcards are used.")
    print("One or more output file names may be specified (the same number")
    print("  as the input file names).")


def x1d(input, output="",
        backcorr="perform", ctecorr="perform", dispcorr="perform",
        helcorr="perform", fluxcorr="perform",
        sporder=None, a2center=None, maxsrch=None,
        globalx=False, extrsize=None,
        bk1size=None, bk2size=None, bk1offst=None, bk2offst=None, bktilt=None,
        backord=None, bksmode="median", bksorder=3,
        blazeshift=None, algorithm="unweighted", xoffset=None,
        verbose=False, timestamps=False, trailer="",
        print_version=False, print_revision=False):
    """Extract a 1-D spectrum from an flt or crj file.

    Parameters
    ----------
    input: str
        Name of the input raw file.

    output: str
        Name of the output file, or "" (the default).  If no name was
        specified, the output name will be constructed from the input name.

    backcorr: str
        If "perform", subtract the background.

    ctecorr: str
        If "perform", apply CTE correction (CCD only).

    dispcorr: str
        If "perform", compute wavelengths from the dispersion relation.

    helcorr: str
        If "perform", correct for heliocentric Doppler shift.

    fluxcorr: str
        If "perform", convert to absolute flux.

    sporder: int or None
        The number of the spectral order to extract.

    a2center: float or None

    maxsrch: float or None

    globalx: bool
        If True, use the global cross correlation offset (i.e. average for
        all orders) for all spectral orders.

    extrsize:  float or None
        Size of extraction box.  None means extrsize is not specified.

    bk1size:  float or None
        Size of first background region.  None means bk1size is not specified.

    bk2size:  float or None
        Size of second background region.  None means bk2size is not specified.

    bk1offst:  float or None
        Offset of first background region.  None means bk1offst is not
                specified.

    bk2offst:  float or None
        Offset of first background region.  None means bk2offst is not
        specified.

    bktilt:  float or None
        Background tilt.  None means bktilt is not specified.

    backord:  int or None
        Background order (0 or 1).  None means backord is not specified.

    bksmode: str
        Background smoothing mode ("off", "median" (the default), or
        "average").

    bksorder: int
        Background smoothing polynomial order (default is 3).

    blazeshift: float or None
        Blaze shift (in pixels).  None means blazeshift is not specified.

    algorithm: str
        Extraction algorithm ("unweighted" (the default) or "sc2d")

    xoffset: float
        Offset in X for slitless extraction.

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
        1 is returned if cs6.e (the calstis host executable) returned a
        non-zero status.  If verbose is True, the value returned by cs6.e
        will be printed.
        2 is returned if the specified input file or files were not found,
        or the numbers of input and output files (if the latter was
        specified) are not the same.
    """

    if print_version:
        status = subprocess.call(["cs6.e", "--version"])
        return 0
    if print_revision:
        status = subprocess.call(["cs6.e", "-r"])
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

        arglist = ["cs6.e"]

        arglist.append(infile)
        if outfiles:
            arglist.append(outfiles[i])

        if verbose:
            arglist.append("-v")
        if timestamps:
            arglist.append("-t")
        if globalx:
            arglist.append("-g")

        switch_was_set = False
        if backcorr == "perform":
            arglist.append("-back")
            switch_was_set = True
        if ctecorr == "perform":
            arglist.append("-cte")
            switch_was_set = True
        if dispcorr == "perform":
            arglist.append("-disp")
            switch_was_set = True
        if helcorr == "perform":
            arglist.append("-hel")
            switch_was_set = True
        if fluxcorr == "perform":
            arglist.append("-flux")
            switch_was_set = True
        if not switch_was_set:
            arglist.append("-x1d")

        if sporder is not None:
            arglist.append("-s")
            arglist.append("%d" % sporder)

        if a2center is not None:
            arglist.append("-c")
            arglist.append("%.10g" % a2center)

        if maxsrch is not None:
            arglist.append("-r")
            arglist.append("%.10g" % maxsrch)

        if extrsize is not None:
            arglist.append("-x")
            arglist.append("%.10g" % extrsize)

        if bk1size is not None:
            arglist.append("-b1")
            arglist.append("%.10g" % bk1size)

        if bk2size is not None:
            arglist.append("-b2")
            arglist.append("%.10g" % bk2size)

        if bk1offst is not None:
            arglist.append("-o1")
            arglist.append("%.10g" % bk1offst)

        if bk2offst is not None:
            arglist.append("-o2")
            arglist.append("%.10g" % bk2offst)

        if bktilt is not None:
            arglist.append("-k")
            arglist.append("%.10g" % bktilt)

        if backord is not None:
            arglist.append("-n")
            arglist.append("%d" % backord)

        if blazeshift is not None:
            arglist.append("-bs")
            arglist.append("%.10g" % blazeshift)

        if bksmode:
            if bksmode == "off":
                arglist.append("-bn")
            elif bksmode == "median":
                arglist.append("-bm")
                arglist.append("-bo")
                arglist.append("%d" % bksorder)
            elif bksmode == "average":
                arglist.append("-bb")
                arglist.append("-bo")
                arglist.append("%d" % bksorder)
            else:
                raise RuntimeError("bksmode must be one of 'off',"
                                   " 'median', 'average'; you specified '%s'" % bksmode)

        if algorithm:
            if algorithm == "unweighted":
                arglist.append("-a")
                arglist.append("unweighted")
            elif algorithm == "sc2d":
                arglist.append("-a")
                arglist.append("unweighted")
                arglist.append("-idt")
            else:
                raise RuntimeError("algorithm must be either 'unweighted'"
                                   " or 'sc2d'; you specified '%s'" % algorithm)

        if xoffset is not None:
            arglist.append("-st")
            arglist.append("%.10g" % xoffset)

        if verbose:
            print("Running x1d on {}".format(infile))
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
    """Return documentation on the x1d function."""
    return x1d.__doc__


def run(configobj=None):
    """TEAL interface for the x1d function."""
    x1d(input=configobj["input"],
        output=configobj["output"],
        backcorr=configobj["backcorr"],
        ctecorr=configobj["ctecorr"],
        dispcorr=configobj["dispcorr"],
        helcorr=configobj["helcorr"],
        fluxcorr=configobj["fluxcorr"],
        sporder=configobj["sporder"],
        a2center=configobj["a2center"],
        maxsrch=configobj["maxsrch"],
        globalx=configobj["globalx"],
        extrsize=configobj["extrsize"],
        bk1size=configobj["bk1size"],
        bk2size=configobj["bk2size"],
        bk1offst=configobj["bk1offst"],
        bk2offst=configobj["bk2offst"],
        bktilt=configobj["bktilt"],
        backord=configobj["backord"],
        bksmode=configobj["bksmode"],
        bksorder=configobj["bksorder"],
        blazeshift=configobj["blazeshift"],
        algorithm=configobj["algorithm"],
        xoffset=configobj["xoffset"],
        verbose=configobj["verbose"],
        timestamps=configobj["timestamps"],
        trailer=configobj["trailer"],
        print_version=configobj["print_version"],
        print_revision=configobj["print_revision"])

if __name__ == "__main__":

    main(sys.argv[1:])
