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
        (options, pargs) = getopt.getopt(args, "rtv:",
                                         ["version"])
    except Exception, error:
        prtOptions()
        sys.exit()

    output = None
    outblev = None
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

def basic2d(input, output=None, outblev=None,
           verbose=False, timestamps=False,
           dqicorr=True, blevcorr=True, doppcorr=True,
           lorscorr=True, glincorr=True, lflgcorr=True,
           biascorr=True, darkcorr=True, flatcorr=True,
           photcorr=True, statflag=True,
           darkscale="",
           trailer=None, print_version=False, print_revision=False):

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

    if output is None:
        outfiles = None
    else:
        outfiles = []
        output1 = output.split()
        for out1 in output1:
            output2 = out1.split(",")
            for out2 in output2:
                outfiles.append(out2)

    if outblev is None:
        outblev_txt = None
    else:
        outblev_txt = []
        outblev1 = outblev.split()
        for out1 in outblev1:
            outblev2 = out1.split(",")
            for out2 in outblev2:
                outblev_txt.append(out2)

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

if __name__ == "__main__":

    main(sys.argv[1:])
