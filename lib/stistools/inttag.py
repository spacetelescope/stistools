#! /usr/bin/env python
import numpy as np
from astropy.io import fits
import astropy.stats
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
    cenx = tag_hdr[0].header['CENTERA1']  # xcenter in c code
    ceny = tag_hdr[0].header['CENTERA2']  # ycenter in c code
    siz_axx = tag_hdr[0].header['SIZAXIS1']  # nx in c code
    siz_axy = tag_hdr[0].header['SIZAXIS2']  # ny in c code
    tzero_mjd = tag_hdr[1].header['EXPSTART']  # MJD zero point

    xcorner = ((cenx - siz_axx / 2.) - 1) * 2
    ycorner = ((ceny - siz_axy / 2.) - 1) * 2

    # Adjust axis sizes for highres, determine binning
    bin_n = 2
    if highres:
        siz_axx *= 2
        siz_axy *= 2
        bin_n = 1

    ltvx = ((bin_n - 2.) / 2. - xcorner) / bin_n
    ltvy = ((bin_n - 2.) / 2. - ycorner) / bin_n
    ltm = 2. / bin_n

    # Read in start and stop time parameters
    if starttime is None:
        starttime = gti_start  # The first START time in the GTI (or first event)

    if increment is None:
        increment = (gti_stop - gti_start)/rcount
    stoptime = starttime + increment

    imset_hdr_ver = 0  # output header value corresponding to imset
    texptime = 0  # total exposure time
    hdu_list = []
    for imset in range(rcount):

        # Get Exposure Times
        #exp_time, expstart, expstop = exp_range(starttime, stoptime, gti_data, tzero_mjd)
        exp_time, expstart, expstop, good_events = exp_range_py(starttime, stoptime, events_data, gti_data, tzero_mjd)

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
        texpend = expstop  # texpend will be expstop of last imset

        if verbose:
            print("imset: {}, start: {}, stop: {}, exposure time: {}".format(imset_hdr_ver,
                                                                            starttime,
                                                                            stoptime,
                                                                            exp_time))

        # Convert events table to accum image
        accum = events_to_accum(good_events, siz_axx, siz_axy, highres)

        # Calculate errors from accum image
        # Note: C takes the square root of the counts, inttag.py uses a more robust confidence interval
        conf_int = astropy.stats.poisson_conf_interval(accum, interval='sherpagehrels', sigma=1)
        err = conf_int[1] - accum  # Error is the difference between upper confidence boundary and the data

        # Copy EVENTS extension header to SCI, ERR, DQ extensions
        sci_hdu = fits.ImageHDU(data=accum, header=tag_hdr[1].header.copy(), name='SCI')
        err_hdu = fits.ImageHDU(data=err, header=tag_hdr[1].header.copy(), name='ERR')
        dq_hdu = fits.ImageHDU(header=tag_hdr[1].header.copy(), name='DQ')

        # Populate extensions
        for hdu in [sci_hdu, err_hdu, dq_hdu]:
            hdu.header['EXPTIME'] = exp_time
            hdu.header['EXPSTART'] = expstart
            hdu.header['EXPEND'] = expstop
            hdu.header['EXTVER'] = imset_hdr_ver

            # Check if image-specific WCS keywords already exist in the tag file (older tag files do)
            keyword_list = list(hdu.header.keys())
            if not any("CTYPE" in keyword for keyword in keyword_list):
                n, k = [keyword[-1] for keyword in keyword_list if "TCTYP" in keyword]
                # Rename keywords
                for val, i in zip([n, k], ['1', '2']):
                    hdu.header.rename_keyword('TCTYP' + val, 'CTYPE' + i)
                    hdu.header.rename_keyword('TCRPX' + val, 'CRPIX' + i)
                    hdu.header.rename_keyword('TCRVL' + val, 'CRVAL' + i)
                    hdu.header.rename_keyword('TCUNI' + val, 'CUNIT' + i)
                hdu.header.rename_keyword('TC{}_{}'.format(n, n), 'CD{}_{}'.format(1, 1))
                hdu.header.rename_keyword('TC{}_{}'.format(n, k), 'CD{}_{}'.format(1, 2))
                hdu.header.rename_keyword('TC{}_{}'.format(k, n), 'CD{}_{}'.format(2, 1))
                hdu.header.rename_keyword('TC{}_{}'.format(k, k), 'CD{}_{}'.format(2, 2))

            # Time tag events table keywords
            hdu.header['WCSAXES'] = 2
            hdu.header['LTM1_1'] = ltm
            hdu.header['LTM2_2'] = ltm
            hdu.header['LTV1'] = ltvx
            hdu.header['LTV2'] = ltvy

            if not highres:
                pass
                hdu.header['CD1_1'] *= 2
                hdu.header['CD1_2'] *= 2
                hdu.header['CD2_1'] *= 2
                hdu.header['CD2_2'] *= 2
                hdu.header['CRPIX1'] = (hdu.header['CRPIX1'] + 0.5) / 2.
                hdu.header['CRPIX2'] = (hdu.header['CRPIX2'] + 0.5) / 2.

        hdu_list.append(sci_hdu)
        hdu_list.append(err_hdu)
        hdu_list.append(dq_hdu)

        starttime = stoptime
        stoptime += increment
        texptime += exp_time

    pri_hdu = fits.PrimaryHDU(header=tag_hdr[0].header.copy())
    pri_hdu.header['NEXTEND'] = imset_hdr_ver * 3  # Three extensions per imset (SCI, ERR, DQ)
    pri_hdu.header['NRPTEXP'] = imset_hdr_ver
    pri_hdu.header['TEXPSTRT'] = texpstart
    pri_hdu.header['TEXPEND'] = texpend
    pri_hdu.header['TEXPTIME'] = texptime

    # Write output file
    hdu_list = [pri_hdu] + hdu_list
    out_hdul = fits.HDUList(hdu_list)
    out_hdul.writeto(output, overwrite=True)


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

    good_events = imset_events[gti_mask]  # All events in the imset within the GTI(s)
    expstart = tzero_mjd + good_events['TIME'][0]/sec_per_day  # exposure start in MJD for imset
    expstop = tzero_mjd + good_events['TIME'][-1] / sec_per_day  # exposure stop in MJD for imset
    exp_time = imset_events['TIME'][-1] - imset_events['TIME'][0] - masked_time  # exposure time in seconds

    return exp_time, expstart, expstop, good_events


def events_to_accum(events_table, size_x, size_y, highres):
    print(size_x, size_y)
    axis1 = events_table['AXIS1']
    axis2 = events_table['AXIS2']
    if highres:
        range_y = size_y
        range_x = size_x
    else:
        range_y = size_y * 2
        range_x = size_x * 2
    accum, xedges, yedges = np.histogram2d(axis2, axis1, bins=[size_y, size_x], range=[[1, range_y], [1, range_x]])
    return accum


if __name__ == "__main__":
    inttag("ob3001xqq_tag.fits", "test.fits", rcount=3, verbose=True, highres=False)
