#! /usr/bin/env python


import os
import sys
import getopt
import glob
import subprocess

from stsci.tools import parseinput, teal

__doc__ = """
Add STIS exposures, rejecting cosmic rays.

Examples
--------

In Python without TEAL:

>>> import stistools
>>> stistools.ocrreject.ocrreject("o3tt02020_flt.fits",
...     "o3tt02020_crj.fits", verbose=True, trailer="o3tt02020.trl")

In Python with TEAL:

>>> from stistools import ocrreject
>>> from stsci.tools import teal
>>> teal.teal("ocrreject")

From command line::

% ./ocrreject.py -v -s o3tt02020_flt.fits o3tt02020_crj.fits
% ./ocrreject.py -r
"""

__taskname__ = "ocrreject"
__version__ = "3.4"
__vdate__ = "13-November-2013"
__author__ = "Phil Hodge, STScI, November 2013."


def main(args):

    if len(args) < 2:
        prtOptions()
        print("At least input and output file names must be specified.")
        sys.exit()

    try:
        (options, pargs) = getopt.getopt(args, "rtv:",
                                         ["version"])
    except Exception as error:
        prtOptions()
        sys.exit()

    verbose = False
    timestamps = False

    for i in range(len(options)):
        if options[i][0] == "--version":
            status = subprocess.call(["cs2.e", "--version"])
            return 0
        if options[i][0] == "-r":
            status = subprocess.call(["cs2.e", "-r"])
            return 0
        elif options[i][0] == "-v":
            verbose = True
        elif options[i][0] == "-t":
            timestamps = True

    nargs = len(pargs)
    if nargs != 2:
        prtOptions()
        sys.exit()
    input = pargs[0]
    output = pargs[1]

    status = ocrreject(input, output,
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
    print("Following the options, list one or more input files")
    print("  (enclosed in quotes if more than one file name is specified")
    print("  and/or if wildcards are used) and one output file name.")


def ocrreject(input, output,
              all=True, crrejtab="", scalense="", initgues="",
              skysub="", crsigmas="",
              crradius=None, crthresh=None, badinpdq=None,
              crmask="",
              verbose=False, timestamps=False,
              trailer="", print_version=False, print_revision=False):
    """Find and reject cosmic rays in STIS data.

    Parameters
    ----------
    input: str
        Name of the input file or files.

    output: str
        Name of the output file.  See all for further information.

    all: bool
        If True (the default), combine all input files into one output
        file.  In this case, output should just be one file name.  If
        False, the number of input and output file names must be the same.

    crrejtab: str
        This argument may be used to override the CRREJTAB value in the
        primary headers of the input files.

    scalense: str
        If specified, this overrides SCALENSE in the CRREJTAB.

    initgues: str
        If specified, this overrides INITGUES in the CRREJTAB.  The
        allowed values are "min" and "med" and "".

    skysub: str
        If specified, this overrides SKYSUB in the CRREJTAB.  The
        allowed values are "none", "mode" and "".

    crsigmas: str
        If specified, this overrides CRSIGMAS in the CRREJTAB.  The
        value should be a comma-separated string of one or more
        integer or float values.  For each such value, calstis will
        perform one cosmic-ray-rejection cycle, with the sigma taken
        from the numerical value that was specified.

    crradius: float or None
        If not None, this overrides CRRADIUS in the CRREJTAB.  This is
        the rejection propagation radius in pixels (e.g. 1.5).  After
        finding an outlier (a cosmic ray hit), adjacent pixels can also
        be flagged and excluded.  Neighboring pixels will be rejected if
        their values are discrepant by more than crthresh * sigmas * noise,
        where noise is based on the noise model (i.e. Poisson noise and
        readout noise).

    crthresh: float or None
        If not None, this overrides CRTHRESH in the CRREJTAB.  This is the
        rejection propagation threshold (e.g. 0.8).  If crthresh = 0 then
        all adjacent pixels (see crradius) will be rejected.

    badinpdq: int or None
        If specified, this overrides BADINPDQ in the CRREJTAB.  This is a
        data quality flag (or bitwise OR of flags) to allow rejection of
        pixels in the input images when forming the "guess" image (the
        image with which to compare the input images when looking for
        outliers).

    crmask: str
        If specified, this overrides CRMASK in the CRREJTAB.  crmask =
        "yes" means that the cosmic rays that are detected should be
        flagged in the DQ (data quality) extensions of the input files.

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
        1 is returned if cs2.e (the calstis host executable) returned a
        non-zero status.  If verbose is True, the value returned by cs2.e
        will be printed.
    """

    if print_version:
        status = subprocess.call(["cs2.e", "--version"])
        return 0
    if print_revision:
        status = subprocess.call(["cs2.e", "-r"])
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

    outfiles = []
    output1 = output.split()
    for out1 in output1:
        if out1:
            output2 = out1.split(",")
            for out2 in output2:
                if out2:
                    outfiles.append(out2)

    n_outfiles = len(outfiles)
    if all:
        if n_outfiles != 1:
            print("You specified {} output files; when all is True,".format(
                n_outfiles))
            print("output must be exactly one file name.")
            return 2
    else:
        n_infiles = len(infiles)
        if n_outfiles != n_infiles:
            print("You specified {} input files but {} output files;".format(
                n_infiles, n_outfiles))
            print("the number of input and output files must be the same.")
            return 2

    if trailer:
        if verbose and os.access(trailer, os.F_OK):
            print("Appending to trailer file {}".format(trailer))
        f_trailer = open(trailer, "a")
        fd_trailer = f_trailer.fileno()
    else:
        f_trailer = None
        fd_trailer = None

    optional_args = []
    if crrejtab:
        optional_args.append("-table")
        optional_args.append(crrejtab)
    if scalense:
        optional_args.append("-scale")
        optional_args.append(scalense)
    if initgues:
        optional_args.append("-init")
        optional_args.append(initgues)
    if skysub:
        optional_args.append("-sky")
        optional_args.append(skysub)
    if crsigmas:
        optional_args.append("-sigmas")
        optional_args.append(crsigmas)
    if crradius:
        optional_args.append("-radius")
        optional_args.append("%.10g" % crradius)
    if crthresh:
        optional_args.append("-thresh")
        optional_args.append("%.10g" % crthresh)
    if badinpdq:
        optional_args.append("-pdq")
        optional_args.append("%d" % badinpdq)
    if crmask:
        if crmask == "yes":
            optional_args.append("-crmask")
            optional_args.append("yes")
        elif crmask == "no":
            optional_args.append("-crmask")
            optional_args.append("no")
        else:
            raise RuntimeError("crmask = %s, must be yes or no." % crmask)

    if all:
        arglist = ["cs2.e"]

        infilestr = "%s" % infiles[0]
        n_infiles = len(infiles)
        for i in range(1, n_infiles):
            infilestr += " %s" % infiles[i]
        arglist.append(infilestr)
        arglist.append(output)

        if verbose:
            arglist.append("-v")
        if timestamps:
            arglist.append("-t")
        arglist.extend(optional_args)

        if verbose:
            print("'{}'".format(str(arglist)))
            print("Running ocrreject on {}".format(infilestr))
        del infilestr
        status = subprocess.call(arglist, stdout=fd_trailer,
                                 stderr=subprocess.STDOUT)
        if status and verbose:
            print("Warning:  status = {}".format(status))
        cumulative_status = status

    else:
        cumulative_status = 0
        for (i, infile) in enumerate(infiles):
            arglist = ["cs2.e"]
            arglist.append(infile)
            arglist.append(outfiles[i])

            if verbose:
                arglist.append("-v")
            if timestamps:
                arglist.append("-t")
            arglist.extend(optional_args)

        if verbose:
            print("Running ocrreject on {}".format(infile))
            print("  {}".format(str(arglist)))
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
    """Return documentation on the ocrreject function."""
    return ocrreject.__doc__


def run(configobj=None):
    """TEAL interface for the ocrreject function."""
    ocrreject(input=configobj["input"],
              output=configobj["output"],
              all=configobj["all"],
              crrejtab=configobj["crrejtab"],
              scalense=configobj["scalense"],
              initgues=configobj["initgues"],
              skysub=configobj["skysub"],
              crsigmas=configobj["crsigmas"],
              crradius=configobj["crradius"],
              crthresh=configobj["crthresh"],
              badinpdq=configobj["badinpdq"],
              crmask=configobj["crmask"],
              verbose=configobj["verbose"],
              timestamps=configobj["timestamps"],
              trailer=configobj["trailer"],
              print_version=configobj["print_version"],
              print_revision=configobj["print_revision"])

if __name__ == "__main__":

    main(sys.argv[1:])
