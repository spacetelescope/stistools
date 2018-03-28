#! /usr/bin/env python
import numpy as np
from astropy.io import fits


def inttag(input, output, starttime=None, increment=None,
            rcount=1, highres=False, allevents=False, verbose=True):
    """Convert an events table of TIMETAG into an integrated ACCUM image.

        Parameters
        ----------
        input: str
            nput file that contains TIMETAG event stream. This is ordinarily a
            FITS file containing two tables. The TIMETAG data are in the table
            with EXTNAME = "EVENTS", and the "good time intervals" are in the
            table with EXTNAME = "GTI". If the GTI table is missing or empty,
            all times will be considered "good".

        output: str
            Name of the output FITS file.

        starttime: float
            Start time for integrating events, in units of seconds since the
            beginning of the exposure. The default value of None means that
            the start time will be set to the first START time in the GTI table.

        increment: float
            Time interval in seconds. The default value of None means integrate
            to the last STOP time in the GTI table.

        rcount: int
            Repeat count, the number of output image sets to create. If rcount is
            greater than 1 then increment must also be specified (but starttime may
            still be set to None).

        highres: bool
            Create a high resolution output image? Default is False.

        allevents: bool
            If allevents is set to yes, all events in the input EVENTS table will
            be accumulated into the output image. The TIME column in the EVENTS
            table will only be used to determine the exposure time, and the GTI
            table will be ignored.

        verbose: bool
            Print additional info?

        Returns
        -------

"""
    # Open Input File (_tag.fits)
    tag_hdr = fits.open(input)


    # Read in GTI Table
    gti_data = tag_hdr['GTI'].data
    # C code checked if it could retrieve a GTI table,
    # If not, it used the events table. Not sure if we
    # want to allow this, leaving it out for now

    # Get Header Info
    cen1 = tag_hdr[0].header['CENTERA1']  # xcenter in c code
    cen2 = tag_hdr[0].header['CENTERA2']  # ycenter in c code
    siz_ax1 = tag_hdr[0].header['SIZAXIS1']  # nx in c code
    siz_ax2 = tag_hdr[0].header['SIZAXIS2']  # ny in c code
    tzero_mjd = tag_hdr[1].header['EXPSTART'] # MJD zeropoint

    # Adjust axis sizes for highres
    if highres:
        siz_ax1 *= 2
        siz_ax2 *= 2

    # Read in start and stop time parameters
    if starttime is None:
        starttime = gti_data['START'][0]  # The first START time in the GTI

    if increment is None:
        endtime = gti_data['STOP'][-1]  # The last STOP time in the GTI
    else:
        endtime = starttime + increment

    imset_hdr_ver = 0  # output header value corresponding to imset
    for imset in range(rcount):

        # Get Exposure Times
        exp_time, expstart, expend = exp_range(starttime,endtime,gti_data,tzero_mjd)
        if exp_time <= 0.:
            if verbose:
                print("Skipping imset, due to no overlap with GTI\n", starttime, endtime)
            starttime = endtime
            endtime += increment
            continue

        imset_hdr_ver += 1

        if imset_hdr_ver == 1: # If first science header, texpstart keyword value is expstart
            texpstart = expstart
        texpend = expend # texpend will be expend of last imset

        if verbose:
            print("imset: {}, start: {}, end: {}, exposure time: {}".format(imset_hdr_ver,
                                                                            starttime,
                                                                            endtime,
                                                                            exp_time))





def exp_range(starttime, endtime, gti_data, tzero_mjd):
    """Calculate exposure time, expstart, and expend """
    start_times = gti_data['START']
    stop_times = gti_data['STOP']
    sec_per_day = 86400

    # Figure out start and stop times from GTI(s)
    if (endtime < start_times[0]) or (starttime > stop_times[-1]):
        exp_time = 0
        expstart = tzero_mjd
        expend = tzero_mjd
        return exp_time, expstart, expend

    exp_time = 0
    for i in range(len(gti_data)):
        # For each GTI
        if (endtime < start_times[i]) or (starttime > stop_times[i]):
            # If the imset range ends before the GTI or starts after the GTI,
           continue # Move on
        start_expt = max([start_times[i], starttime])
        end_expt = min([stop_times[i], endtime])
        exp_time += (end_expt - start_expt)
    if exp_time <= 0.0:
        expstart = tzero_mjd
        expend = tzero_mjd
        return exp_time, expstart, expend
    # The following assumes the GTI are sorted by time
    for i in range(len(gti_data)):
        if (starttime >= start_times[i]) and (starttime <= stop_times[i]):
            start_expt = starttime
        elif starttime <= start_times[i]:
            start_expt = start_times[i]
    expstart = tzero_mjd + start_expt/sec_per_day

    for i in range(len(gti_data)):
        if (endtime >= start_times[i]) and endtime <=stop_times[i]:
            end_expt = endtime
        elif endtime < start_times[i]:
            end_expt = stop_times[i]
    expend = tzero_mjd + end_expt/sec_per_day
    return exp_time, expstart, expend












if __name__ == "__main__":
    inttag("o5d601050_tag.fits", "test.fits")
