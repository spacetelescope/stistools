#! /usr/bin/env python
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


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
            to the last STOP time in the GTI table, divided by rcount.

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

    # Read in GTI (Good Time Interval) Table
    gti_data = tag_hdr['GTI'].data
    events_data = tag_hdr[1].data

    # Determine start and stop times for counting events

    # If the user sets the allevents flag, the time interval is determined by the
    # first and last event in the events table. Otherwise, it is determined by the
    # first START and the last STOP in the GTI Table
    if allevents:
        gti_start = tag_hdr[1].data['TIME'][0]
        gti_stop = tag_hdr[1].data['TIME'][-1]
    else:
        gti_start = gti_data['START'][0]
        gti_stop = gti_data['STOP'][-1]

    # C code checked if it could retrieve a GTI table,
    # If not, it used the events table. Not sure if we
    # want to allow this, leaving it out for now

    # Get Header Info
    cen1 = tag_hdr[0].header['CENTERA1']  # xcenter in c code
    cen2 = tag_hdr[0].header['CENTERA2']  # ycenter in c code
    siz_ax1 = tag_hdr[0].header['SIZAXIS1']  # nx in c code
    siz_ax2 = tag_hdr[0].header['SIZAXIS2']  # ny in c code
    tzero_mjd = tag_hdr[1].header['EXPSTART']  # MJD zero point

    # Adjust axis sizes for highres
    if highres:
        siz_ax1 *= 2
        siz_ax2 *= 2

    # Read in start and stop time parameters
    if starttime is None:
        starttime = gti_start  # The first START time in the GTI (or first event)

    if increment is None:
        increment = (gti_stop - gti_start)/rcount
    stoptime = starttime + increment

    imset_hdr_ver = 0  # output header value corresponding to imset
    for imset in range(rcount):

        # Get Exposure Times
        #exp_time, expstart, expstop = exp_range(starttime, stoptime, gti_data, tzero_mjd)
        exp_time, expstart, expstop, good_events = exp_range_py(starttime, stoptime, events_data, gti_data, tzero_mjd)

        #print(exp_time, expstart, expstop, good_events)
        #if exp_time <= 0.:
        if len(good_events) == 0:
            if verbose:
                print("Skipping imset, due to no overlap with GTI\n", starttime, stoptime)
            starttime = stoptime
            stoptime += increment
            continue

        imset_hdr_ver += 1

        if imset_hdr_ver == 1:  # If first science header, texpstart keyword value is expstart
            texpstart = expstart
        texpstop = expstop  # texpstop will be expstop of last imset

        if verbose:
            print("imset: {}, start: {}, stop: {}, exposure time: {}".format(imset_hdr_ver,
                                                                            starttime,
                                                                            stoptime,
                                                                            exp_time))

        accum = events_to_accum(good_events, siz_ax1, siz_ax2)
        plt.imshow(accum,vmin=0,vmax=10,origin='lower')
        plt.show()
        starttime = stoptime
        stoptime += increment


def exp_range(starttime, stoptime, gti_data, tzero_mjd):
    """Calculate exposure time, expstart, and expstop """
    start_times = gti_data['START']
    stop_times = gti_data['STOP']
    sec_per_day = 86400

    # Figure out start and stop times from GTI(s)
    if (stoptime < start_times[0]) or (starttime > stop_times[-1]):  # Outside total GTI Range
        exp_time = 0
        expstart = tzero_mjd
        expstop = tzero_mjd
        return exp_time, expstart, expstop

    exp_time = 0
    for i in range(len(gti_data)):  # For each GTI
        if (stoptime < start_times[i]) or (starttime > stop_times[i]):  # Not in GTI
            continue
        start_expt = max([start_times[i], starttime])
        stop_expt = min([stop_times[i], stoptime])
        exp_time += (stop_expt - start_expt)
    if exp_time <= 0.0:
        expstart = tzero_mjd
        expstop = tzero_mjd
        return exp_time, expstart, expstop

    # The following assumes the GTI are sorted by time
    for i in range(len(gti_data)):  # For each GTI
        if (starttime >= start_times[i]) and (starttime <= stop_times[i]):  # Starts within GTI
            start_expt = starttime
        elif starttime <= start_times[i]:  # Starts before GTI
            start_expt = start_times[i]
    expstart = tzero_mjd + start_expt/sec_per_day

    for i in range(len(gti_data)):  # For each GTI
        if (stoptime >= start_times[i]) and stoptime <= stop_times[i]:  # Stops within GTI
            stop_expt = stoptime
        elif stoptime > stop_times[i]:  # Stops after GTI, different from C code, write a test for this
            stop_expt = stop_times[i]
    expstop = tzero_mjd + stop_expt/sec_per_day
    return exp_time, expstart, expstop


def exp_range_py(starttime, stoptime, events_data, gti_data, tzero_mjd):
    """Calculate exposure time, expstart, and expstop and mask imset"""
    sec_per_day = 86400
    imset_events = events_data[(events_data['TIME'] > starttime) * (events_data['TIME'] < stoptime)]  # within imset
    if len(imset_events) == 0:  # No events in imset
        exp_time = 0
        expstart = tzero_mjd
        expstop = tzero_mjd
        return exp_time, expstart, expstop, imset_events

    # Mask events in imset if there are any lapses in GTI
    gti_mask = [True] * len(imset_events)
    masked_time = 0
    for gti in gti_data:
        mask = (imset_events['TIME'] > gti[0]) * (imset_events['TIME'] < gti[1])
        gti_mask *= mask
        masked_events = imset_events[~gti_mask]

        # Track exposure time lost due to non-GTI
        if len(masked_events) > 0:
            masked_time += masked_events['TIME'][-1] - masked_events['TIME'][0]

    good_events = imset_events[gti_mask] # All events in the imset within the GTI(s)
    expstart = tzero_mjd + good_events['TIME'][0]/sec_per_day  # exposure start in MJD for imset
    expstop = tzero_mjd + good_events['TIME'][-1] / sec_per_day  # exposure stop in MJD for imset
    exp_time = good_events['TIME'][-1] - good_events['TIME'][0] - masked_time  # exposure time in seconds

    return exp_time, expstart, expstop, good_events


def events_to_accum(events_table, size_x, size_y):
    print(size_x, size_y)
    axis1 = events_table['AXIS1']
    axis2 = events_table['AXIS2']
    accum, xedges, yedges = np.histogram2d(axis2, axis1, bins=[size_y, size_x], range=[[0, size_y*2], [0, size_x*2]])
    return accum



if __name__ == "__main__":
    inttag("ob3001xqq_tag.fits", "test.fits", rcount=10, verbose=True, highres = False)
