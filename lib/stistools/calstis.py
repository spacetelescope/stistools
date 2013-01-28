#! /usr/bin/env python

"""
Calibrate STIS data.

Examples
--------

In Python without TEAL:

>>> import stistools
>>> stistools.calstis.calstis("o66p01020_raw.fits", verbose=True,
...                           trailer="o66p01020.trl")

In Python with TEAL:

>>> from stistools import calstis
>>> from stsci.tools import teal
>>> teal.teal("calstis")

In Pyraf:

>>> import stistools
>>> teal calstis

From command line::

% ./calstis.py -v -s o66p01020_raw.fits out/
% ./calstis.py -r

"""

from __future__ import division         # confidence unknown
import os
import sys
import getopt
import glob
import subprocess

from stsci.tools import parseinput,teal

__taskname__ = "calstis"
__version__ = "3.0"
__vdate__ = "14-January-2013"
__author__ = "Phil Hodge, STScI, January 2013."

def main(args):

    if len(args) < 1:
        prtOptions()
        print("At least a raw file name must be specified.")
        sys.exit()

    try:
        (options, pargs) = getopt.getopt(args, "srtvw:",
                                         ["version"])
    except Exception, error:
        prtOptions()
        sys.exit()

    outroot = None
    wavecal = None
    verbose = False
    timestamps = False
    savetmp = False

    for i in range(len(options)):
        if options[i][0] == "--version":
            status = subprocess.call(["cs0.e", "--version"])
            return 0
        if options[i][0] == "-r":
            status = subprocess.call(["cs0.e", "-r"])
            return 0
        elif options[i][0] == "-v":
            verbose = True
        elif options[i][0] == "-t":
            timestamps = True
        elif options[i][0] == "-s":
            savetmp = True
        elif options[i][0] == "-w":
            wavecal = options[i][1]

    nargs = len(pargs)
    if nargs < 1 or nargs > 2:
        prtOptions()
        sys.exit()
    input = pargs[0]
    if nargs == 2:
        outroot = pargs[1]

    status = calstis(input, wavecal=wavecal, outroot=outroot,
                     savetmp=savetmp,
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
    print("  -w wavecal")
    print("")
    print("Following the options, list one or more input raw file names,")
    print("  enclosed in quotes if more than one file name is specified")
    print("  and/or if wildcards are used.")
    print("An output directory (include a trailing '/') or a root name for")
    print("  the output files may be specified.")

def calstis(input, wavecal=None, outroot=None, savetmp=False,
            verbose=False, timestamps=False,
            trailer=None, print_version=False, print_revision=False):
    """Calibrate STIS data.

    Parameters
    ----------
    input: str
        Name of the input raw file.

    wavecal: str or None
        Name of the input wavecal file, or None (the default).  This is
        only needed if the name is not the "normal" name
        (rootname_wav.fits).

    outroot: str or None
        Root name for the output files, or None (the default).  This can
        be a directory name, in which case the string must end in '/'.

    savetmp: bool
        True if calstis should not delete temporary files.

    verbose: bool
        If True, calstis will print more info.

    timestamps: bool
        If True, calstis will print the date and time at various points
        during processing.

    trailer: str or None
        If specified (i.e. if not None), the standard output and standard
        error will be written to this file instead of to the terminal.
        Note, however, that if print_version or print_revision is
        specified, the value will be printed to the terminal, and any
        name given for the trailer will be ignored.

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
        1 is returned if cs0.e (the calstis host executable) returned a
        non-zero status.  If verbose is True, the value returned by cs0.e
        will be printed.
        2 is returned if the specified input file or files were not found.
    """

    if print_version:
        status = subprocess.call(["cs0.e", "--version"])
        return 0
    if print_revision:
        status = subprocess.call(["cs0.e", "-r"])
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

    if trailer is None:
        f_trailer = None
    else:
        f_trailer = open(trailer, "w")

    for infile in infiles:

        arglist = ["cs0.e"]
        if verbose:
            arglist.append("-v")
        if timestamps:
            arglist.append("-t")
        if savetmp:
            arglist.append("-s")
        arglist.append(infile)
        if outroot is not None:
            arglist.append(outroot)
        if wavecal is not None:
            arglist.append("-w")
            arglist.append("%s" % wavecal)

        if verbose:
            print("Running calstis on %s" % infile)
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
    """Return documentation on the calstis function."""
    return calstis.__doc__

def run(configobj=None):
    """TEAL interface for the calstis function."""
    calstis(configobj["input"],
            configobj["wavecal"],
            configobj["outroot"],
            configobj["savetmp"],
            configobj["verbose"],
            configobj["timestamps"],
            configobj["trailer"],
            configobj["print_version"],
            configobj["print_revision"])

if __name__ == "__main__":

    main(sys.argv[1:])
