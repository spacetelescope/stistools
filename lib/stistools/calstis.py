#! /usr/bin/env python

from __future__ import division         # confidence unknown
import os
import sys
import getopt
import glob
import subprocess

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
            verbose=True, timestamps=False,
            trailer=None, print_version=False, print_revision=False):

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

if __name__ == "__main__":

    main(sys.argv[1:])
