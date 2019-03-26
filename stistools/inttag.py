#! /usr/bin/env python
import numpy as np
from astropy.io import fits
import astropy.stats
from astropy import units as u
from astropy.time import Time
from datetime import datetime as dt

__doc__ = """
The task :func:`inttag` converts an events table of TIME-TAG mode STIS data into a raw, time-integrated ACCUM 
image. By default, :func:`inttag` only integrates over the good time intervals (GTI), though the user can 
choose to integrate over the entire exposure time by setting allevents=True. The output image can be 
calibrated as any other raw image.

The input file for :func:`inttag` is an event stream table of TIME-TAG mode produced by generic conversion. 
The data will be Doppler corrected (as required for medium and high resolution spectroscopic modes). This file 
will consist of a primary header with no data, and two binary table extensions. The primary header is identical 
in structure to the primary header of an ACCUM mode image. The first binary table (EXTNAME=EVENTS) contains a 
list of the events themselves (i.e. science data as an event stream), and the second binary table (EXTNAME=GTI) 
contains a list of good time intervals for the TIMETAG exposure. Columns "TIME", "AXIS1", and "AXIS2" in the 
EVENTS table are read. Columns "START" and "STOP" in the GTI table are read.

The output image is a time integrated (ACCUM mode) image with the same structure as any other STIS MAMA raw 
image  (i.e. primary header followed by a single or series of triplet extensions: SCI, ERR, DQ). The number of 
triplets is determined by the value of rcount. The time interval in the Nth triplet covers from 
(starttime + (N-1)*increment) to (starttime + N*increment). The exposure time in each interval need not be 
identical, because events are included in the image only if they occur during "good time intervals" (as 
determined by the GTI extension table). The keyword OBSMODE in the primary header of the output image will 
still be set to "TIME-TAG".

The output science image is ready to be calibrated (see :func:`calstis`, :func:`crreject`, :func:`basic2d`, 
:func:`x2d`, :func:`x1d`). 

Examples
--------

:func:`inttag` with default values:

>>> import stistools
>>> stistools.inttag.inttag("oddv01050_tag.fits", "oddv01050_raw.fits")

:func:`inttag` with highres output:

>>> import stistools
>>> stistools.inttag.inttag("oddv01050_tag.fits", "oddv01050_raw.fits", highres=True)

:func:`inttag` with multiple output imsets (5 count regions of 200s each):

>>> import stistools
>>> stistools.inttag.inttag("oddv01050_tag.fits", "oddv01050_raw.fits", rcount = 5, increment = 200)

"""

__taskname__ = "inttag"
__version__ = "1.0"
__vdate__ = "13-November-2018"
__author__ = "Python: Doug Branton, C code: R. Katsanis, N. Zarate, Phil Hodge"


def inttag(tagfile, output, starttime=None, increment=None,
           rcount=1, highres=False, allevents=False, verbose=True):
    """Convert an events table of TIMETAG into an integrated ACCUM image.

        Parameters
        ----------
        tagfile: str
            input file that contains TIMETAG event stream. This is ordinarily a
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
            greater than 1 and increment is not specified, will subdivide the total exposure time by rcount.

        highres: bool
            Create a high resolution output image? Default is False.

        allevents: bool
            If allevents is set to True, all events in the input EVENTS table will
            be accumulated into the output image. The TIME column in the EVENTS
            table will only be used to determine the exposure time, and the GTI
            table will be ignored.

        verbose: bool
            Print additional info?

        Returns
        -------

"""
    # Open Input File (_tag.fits)
    with fits.open(tagfile) as tag_hdr:
        # Read in TIME-TAG data (Events data and GTI data)
        events_data = tag_hdr[1].data
        if allevents:  # If allevents, ignore GTI and generate gti_data based on the time of the first and last event
            gti_data = np.rec.array([(events_data['TIME'][0], events_data['TIME'][-1])],
                                    formats=">f8,>f8", names='START, STOP')
        else:  # Otherwise, retrieve the GTIs from the GTI extension
            gti_data = tag_hdr['GTI'].data

        # Read in relevant header info
        tag_pri_hdr = tag_hdr[0].header
        cenx = tag_pri_hdr['CENTERA1']  # xcenter in c code
        ceny = tag_pri_hdr['CENTERA2']  # ycenter in c code
        siz_axx = tag_pri_hdr['SIZAXIS1']  # nx in c code
        siz_axy = tag_pri_hdr['SIZAXIS2']  # ny in c code
        tag_sci_hdr = tag_hdr[1].header
        tzero_mjd = tag_sci_hdr['EXPSTART']  # MJD zero point

    # Determine start and stop times for counting events
    gti_start = gti_data['START'][0]
    gti_stop = gti_data['STOP'][-1]

    # Calculate corners from array size and centers
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
    if starttime is None or starttime < gti_start:
        starttime = gti_start  # The first START time in the GTI (or first event)

    if increment is None:
        increment = (gti_stop - gti_start)/rcount

    stoptime = starttime + increment

    imset_hdr_ver = 0  # output header value corresponding to imset
    texptime = 0  # total exposure time
    hdu_list = []
    for imset in range(rcount):
        # Truncate stoptime at last available event time (GTI or allevents) if it exceeds that
        if stoptime > gti_stop:
            stoptime = gti_stop
        # Get Exposure Times
        exp_time, expstart, expstop, good_events = exp_range(starttime, stoptime, events_data, gti_data, tzero_mjd)
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
            print("imset: {}, start: {}, stop: {}, exposure time: {}".format(imset_hdr_ver, starttime, stoptime,
                                                                             exp_time))

        # Convert events table to accum image
        accum = events_to_accum(good_events, siz_axx, siz_axy, highres)

        # Calculate errors from accum image
        # Note: C version takes the square root of the counts, inttag.py uses a more robust confidence interval
        conf_int = astropy.stats.poisson_conf_interval(accum, interval='sherpagehrels', sigma=1)
        err = conf_int[1] - accum  # error is the difference between upper confidence boundary and the data

        # Copy EVENTS extension header to SCI, ERR, DQ extensions
        sci_hdu = fits.ImageHDU(data=accum, header=tag_sci_hdr.copy(), name='SCI')
        err_hdu = fits.ImageHDU(data=err, header=tag_sci_hdr.copy(), name='ERR')
        dq_hdu = fits.ImageHDU(header=tag_sci_hdr.copy(), name='DQ')

        # Generate datetime for 'DATE' header keyword
        dtstr = str(dt.utcnow())
        date, h, m, s = [dtstr.split()[0], dtstr.split()[1].split(':')[0], dtstr.split()[1].split(':')[1],
                         str(round(float(dtstr.split()[1].split(':')[-1])))]
        if len(s) == 1:
            s = '0' + s

        dtval = date + 'T' + h + ':' + m + ':' + s

        # Populate extensions
        for idx, hdu in enumerate([sci_hdu, err_hdu, dq_hdu]):

            hdu.header['EXPTIME'] = exp_time
            hdu.header['EXPSTART'] = expstart
            hdu.header['EXPEND'] = expstop

            date_obs, time_obs = Time(float(expstart), format='mjd').isot.split('T')
            hdu.header['DATE-OBS'] = date_obs
            hdu.header['TIME-OBS'] = time_obs

            hdu.header['EXTVER'] = imset_hdr_ver
            hdu.header['DATE'] = (dtval, "Date FITS file was generated")
            hdu.header['ORIGIN'] = "stistools inttag.py"

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

            # Convert keyword values to lowres scale if not highres
            if not highres:
                hdu.header['CD1_1'] *= 2
                hdu.header['CD1_2'] *= 2
                hdu.header['CD2_1'] *= 2
                hdu.header['CD2_2'] *= 2
                hdu.header['CRPIX1'] = (hdu.header['CRPIX1'] + 0.5) / 2.
                hdu.header['CRPIX2'] = (hdu.header['CRPIX2'] + 0.5) / 2.

            # Populate DQ header with dq specific keywords
            if idx == 2:
                hdu.header['NPIX1'] = siz_axx
                hdu.header['NPIX2'] = siz_axy
                hdu.header['PIXVALUE'] = 0  # Fixes issue with calstis not running on raw output files

        # Append imset extensions to header list
        hdu_list.append(sci_hdu)
        hdu_list.append(err_hdu)
        hdu_list.append(dq_hdu)

        # Prepare start and stop times for next image in imset
        starttime = stoptime
        stoptime += increment
        texptime += exp_time

    # Copy tag file primary header to output header
    pri_hdu = fits.PrimaryHDU(header=tag_pri_hdr.copy())

    # Add/Modify primary header keywords
    pri_hdu.header['NEXTEND'] = imset_hdr_ver * 3  # Three extensions per imset (SCI, ERR, DQ)
    pri_hdu.header['NRPTEXP'] = imset_hdr_ver
    pri_hdu.header['TEXPSTRT'] = texpstart
    pri_hdu.header['TEXPEND'] = texpend
    pri_hdu.header['TEXPTIME'] = texptime
    pri_hdu.header['BINAXIS1'] = bin_n
    pri_hdu.header['BINAXIS2'] = bin_n
    pri_hdu.header['FILENAME'] = output.split('/')[-1]
    pri_hdu.header['DATE'] = (dtval, "Date FITS file was generated")
    pri_hdu.header['ORIGIN'] = "stistools inttag.py"

    if not highres:
        pri_hdu.header['LORSCORR'] = "COMPLETE"  # Corr flag detailing MAMA data conversion to low res

    # Write output file
    hdu_list = [pri_hdu] + hdu_list
    out_hdul = fits.HDUList(hdu_list)
    out_hdul.writeto(output, overwrite=True)


def exp_range(starttime, stoptime, events_data, gti_data, tzero_mjd):
    """Calculate exposure time, expstart, and expstop and mask imset

        Parameters
        ----------
        starttime: float
            Start time of the imset in seconds

        stoptime: float
            Stop time of the imset in seconds

        events_data: record array
            Record array of timetag events.

        gti_data: record array
            Record array of good time intervals (GTIs).

        tzero_mjd: bool
            Modified Julian Date (MJD) corresponding to the beginning of the exposure

        Returns
        -------
        exp_time: float
            Total exposure time in seconds for the given imset. This number accounts for any exposure time lost
            to non-GTI time (if the user is not using allevents).

        expstart: float

            Start time of the imset exposure

        expstop: float

            Stop time of the imset exposure

        good_events: float
            The events list within the imset exposure time and within the GTIs.


"""
    sec_per_day = (1*u.day).to(u.second).value
    imset_events = events_data[(events_data['TIME'] > starttime) * (events_data['TIME'] < stoptime)]  # within imset

    if len(imset_events) == 0:  # No events in imset
        exp_time = 0
        expstart = tzero_mjd
        expstop = tzero_mjd
        return exp_time, expstart, expstop, imset_events

    # Mask events in imset if there are any lapses in GTI
    gti_mask = np.array([False] * len(imset_events))  # Start by assuming all events are outside all GTIs
    for gti in gti_data:
        # Create mask of events within GTI
        mask = (imset_events['TIME'] > gti[0]) * (imset_events['TIME'] < gti[1])
        gti_mask = np.logical_or(gti_mask, mask)  # OR global gti mask with local gti mask

    good_events = imset_events[gti_mask]  # All events in the imset within the GTI(s)
    if len(good_events) == 0:
        exp_time = 0
        expstart = tzero_mjd
        expstop = tzero_mjd
        return exp_time, expstart, expstop, good_events
    expstart = tzero_mjd + good_events['TIME'][0] / sec_per_day  # exposure start in MJD for imset
    expstop = tzero_mjd + good_events['TIME'][-1] / sec_per_day  # exposure stop in MJD for imset

    # Determine GTI gap regions
    gaps = []
    if len(gti_data) > 1:
        for i, gti in enumerate(gti_data):
            if i == 0:
                continue
            gaps.append((gti_data[i - 1][1], gti[0]))

    # Calculate exposure time lost due to non-GTI overlap
    exptime_loss = 0
    for gap in gaps:
        if gap[1] <= stoptime and gap[0] >= starttime:
            exptime_loss += gap[1] - gap[0]
        elif gap[1] >= stoptime and gap[0] <= stoptime:
            exptime_loss += stoptime - gap[0]
        elif gap[1] >= starttime and gap[0] <= starttime:
            exptime_loss += gap[1] - starttime
        else:
            continue

    exp_time = stoptime - starttime - exptime_loss  # exposure time in seconds

    return exp_time, expstart, expstop, good_events


def events_to_accum(events_data, size_x, size_y, highres):
    """Map timetag events to a 2d accum image array.

        Parameters
        ----------
        events_data: record array
            Record array of timetag events.

        size_x: int
            Number of pixels on axis 1 of the detector.

        size_y: int
            Number of pixels on axis 2 of the detector.

        highres: bool
            Boolean value indicating whether the output accum image is in high or low resolution.

        Returns
        -------
        accum: array
            2d image of all events in the imset on the detector.


"""
    # Extract (x,y) event locations from events_table
    axis1 = events_data['AXIS1']
    axis2 = events_data['AXIS2']

    # Determine resolution-appropriate binning
    if highres:
        range_y = size_y
        range_x = size_x
    else:
        range_y = size_y * 2
        range_x = size_x * 2

    # Map events to an accum image using a 2d histogram
    accum, xedges, yedges = np.histogram2d(axis2, axis1, bins=[size_y, size_x], range=[[1, range_y], [1, range_x]])

    return accum

