#!/usr/bin/env python

import warnings
import time
import datetime
import os
import shutil
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
from astropy.table import Table, vstack
from astropy import units as u, constants as c
from astropy.units.core import UnitsWarning
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
from astropy.coordinates import ITRS, GCRS, ICRS
from astroquery.jplhorizons import Horizons


__doc__ = """
This task :func:`barycentric_correction` calculates the barycentric correction
between HST and the Solar System barycenter for times in HST/STIS observations.
It updates the times in uncalibrated (``tag``, ``raw``) and calibrated data. It
stores the old times in a new column for tag data and in new header keywords
for raw and calibrated data. The function puts the various delay terms and
positions in the header as well.

The function works by first determining the position of the Earth with
Astropy and HST with either JPL Horizons or the STScI-provided HST
orbital files.

For ``*tag.fits`` files with many times to convert, performance can be slow due
to the high number of coordinate transformations.

Stand-alone tasks :func:`calc_delay_jpl` :func:`calc_delay_orbfile` can be used
to calculate barycentric corrections without file modifications.

Task :func:`combine_hst_orb_files` combines two neighboring HST orbital files
together for use on MAST data spanning their coverage.

Task :func:`odelay_file_compare` shows differences in times between two STIS
FITS files. This is useful in what barycentric corrections were calculated.

Examples
--------

:func:`barycentric_correction` with default values:

>>> import stistools
>>> stistools.barycentric_correction.barycentric_correction("oep502040_x1d.fits")
"""

__taskname__ = "barycentric_correction"
__version__ = "1.0"
__vdate__ = "25-March-2026"
__author__ = "J. Lothringer"

SECPERDAY = 86400  # number of sec in a day

warnings.simplefilter('ignore', category=UnitsWarning)


class OrbFileError(ValueError):
    """Orbit file does not cover the requested time range.
    """
    pass


def barycentric_correction(table_names, verbose=True, distance=1e9,
                           hst_orb=None, time_script=False, outfiles=None,
                           time_system='TDB'):
    """
    Calculates time-delay barycentric corrections from HST's position
    to the Solar System barycenter. This correction includes the classic
    geometric Roemer delay, as well as the general relativistic Einstein
    delay.

    This function uses modern Astropy tools to calculate the position
    of the barycenter, Earth's location, and deal with the time standards.
    The function has been tested against a Python implementation of
    the previous `odelaytime` IRAF task and found to be consistent.
    We have also tested against other barycentric correction tools,
    including `barycorr`, `astroutils`, and `pintbary`.

    HST's changing location around the Earth can lead to time-delay
    differences of up to ~46 milliseconds. HST's location can be
    determined either through STScI-provided `HST orbital files <https://www.stsci.edu/~STIS/monitors/ephemeris_files.html>`_,
    or through a query to JPL Horizons.

    These calculations are accurate to within 1 millisecond
    outside the Solar System, and to within 5 milliseconds inside
    the Solar System.

    Parameters
    ----------
    table_names: str or Iterable[str]
        List of strings with the file names to be time-corrected.

    verbose: bool
        Prints completion messages during execution. Default=True

    distance: float
        Distance the object is from HST in AU. Default is a trillion
        AU. Most important for objects in our Solar System as it
        is repsonsible for second-order correction, up to minutes.
        At 1 parsec, the correction can be on the order of a few ms.

    hst_orb: str or None
        Name of HST orbital file (generally starts p, ends as .fit) that
        covers the time of the observations. If not provided, JPL Horizons
        will be used to get HST's orbital position.

    time_script: bool
        Set to True if you want to time how long the script takes. Useful
        for debugging, especially for ``tag`` files with large numbers of
        events.

    outfiles: str or list[str] or None
        If not None, then it is a list of output files for each table_names.
        Each table_name will be copied over to the corresponding outfile name.
        Default None.

    time_system: str
        Define as either ``"TDB"`` or ``"UTC"`` for final time standard conversion.
        They will be different by about 69.184 seconds plus or minus a few
        ms depending on where Earth is in its orbit. Final results will
        then be in either BJD_TDB or BJD_UTC. Default is ``"TDB"``.

    Returns
    -------
    None
        However, the file is written to output.

    Notes
    -----
    If an :class:`OrbFileError` is raised, you may need to combine two
    neighboring HST orbital files using :func:`combine_hst_orb_files`.
    """

    if time_system not in ['UTC', 'TDB']:
        print('time_system not recognized.')
        print('Use "UTC" or "TDB".')
        raise ValueError

    if time_script:
        tstart = time.time()

    # Allow input of one filename as a str:
    if isinstance(table_names, (str, bytes)):
        table_names = [table_names]
    if outfiles is not None:
        if isinstance(outfiles(str, bytes)):
            outfiles = [outfiles]
        if len(table_names) != len(outfiles):
            raise ValueError("If 'outfiles' is specified, it must be the same length as 'table_names'.")

    for ii, in_table_file in enumerate(table_names):

        if verbose:
            print(f"\nbarycentric_correction: processing {in_table_file} ...")

        # Copy to new outfile, otherwise overwrite input file.
        if outfiles is None:
            filename = in_table_file
        else:
            filename = outfiles[ii]
            if verbose:
                print(f"Copying {in_table_file} to {filename}")
            shutil.copy(in_table_file, filename)

        in_hdul = fits.open(filename, mode='update')

        # determine the file type, based on the first extension
        # original code (ln 124) contains some kind of catch all
        # for SCI_IMAGE? It also caught anything outside of those options,
        # might want to add that back in as well.
        extname = in_hdul[1].header['EXTNAME']
        if extname == "EVENTS":
            filetype = "EVENTS_TABLE"
        elif extname == "SCI":
            filetype = "X1D_TABLE"
        else:
            raise ValueError(f"Unexpected extension name:  {extname}")

        if in_hdul[0].header.get("DLAYCORR", "PERFORM") == "COMPLETE":
            print(f"{in_table_file} has already been corrected, so no further processing"
                  " will be applied to this file")
            continue

        # COS has no TEXPSTRT in primary header, so use EXPSTART in first ext
        if in_hdul[0].header['INSTRUME'] == "STIS":
            mjd1 = in_hdul[0].header['TEXPSTRT']
            mjd2 = in_hdul[0].header['TEXPEND']
        elif in_hdul[0].header['INSTRUME'] == "COS":
            mjd1 = in_hdul[1].header['EXPSTART']
            mjd2 = in_hdul[1].header['EXPEND']
        else:
            raise ValueError('baycentric_correction only works with STIS and COS files.')

        ra = in_hdul[0].header['RA_TARG']
        dec = in_hdul[0].header['DEC_TARG']

        if 'NEXTEND' in in_hdul[0].header:
            nextend = in_hdul[0].header['NEXTEND']
        else:
            nextend = 1

        if filetype == "EVENTS_TABLE":
            # assume (for now) that the last extension is a GTI table
            nevents_tab = max(nextend - 1, 1)
            nsci_ext = 0
        elif filetype == "X1D_TABLE":
            nevents_tab = 0
            nsci_ext = nextend
        else:
            nevents_tab = 0
            nsci_ext = nextend

        # Get the delay at the exposure start time.  We'll subtract this
        # delay from each time that we update, if that time is relative
        # to EXPSTART (which must be the same as TEXPSTRT).
        # t0_delay = all_delay(mjd1, parallax, objvec, earth_ephem_table,
        #                     len(earth_ephem_table), obs_ephem_table, npts_obs)

        if hst_orb is None:
            t0_delay = calc_delay_jpl(mjd1, ra, dec, distance=distance)
        else:
            t0_delay = calc_delay_orbfile(mjd1, ra, dec, hst_orb, distance=distance)

        if filetype == "EVENTS_TABLE":
            # update times in GTI table
            if 'GTI' in in_hdul:
                gti_tab = in_hdul['GTI'].data
                for row in gti_tab:
                    for keyword in ['START', 'STOP']:
                        tm = row[keyword]
                        epoch = mjd1 + tm / SECPERDAY

                        if hst_orb is None:
                            delta_sec = calc_delay_jpl(epoch, ra, dec, distance=distance)
                        else:
                            delta_sec = calc_delay_orbfile(epoch, ra, dec, hst_orb, distance=distance)

                        row[keyword] = tm + (delta_sec - t0_delay).value * SECPERDAY

                in_hdul.flush()

                if verbose:
                    print("GTI extension has been updated")

            else:
                # If there's no GTI, last extensions is an EVENTS table
                nevents_tab += 1

        if time_script:
            tcheck1 = time.time()
            print(f'Checkpoint 1: {tcheck1 - tstart} s')

        # loop through all EVENTS extensions (there will be none if
        # the input is an x1d or image file)
        for e_indx in range(1, nevents_tab + 1):
            in_col = [x.name for x in in_hdul['EVENTS', e_indx].data.columns if x.name.lower().strip() == 'time'][0]

            events_tab = in_hdul['EVENTS', e_indx]
            mjd1 = events_tab.header['EXPSTART']
            mjd2 = events_tab.header['EXPEND']

            # find the time (since EXPSTART) column in the table
            nrows = len(events_tab.data[in_col])

            # pull out just the time array
            time_array = events_tab.data[in_col]

            # go through each row
            if verbose:
                print(f"Number of times in extension: {nrows}")

            # re-written to do all rows at once, so we can interpolate times,
            # instead of calling Horizons 4.5 million times
            epoch_array = mjd1 + time_array / SECPERDAY
            if hst_orb is None:
                delta_sec = calc_delay_jpl(epoch_array, ra, dec, distance=distance)
            else:
                delta_sec = calc_delay_orbfile(epoch_array, ra, dec, hst_orb, distance=distance)

            # This is correcting for the difference between the new time interval
            # and the *corrected* exposure start time
            # (i.e., the add'l correction now that the Earth and HST have moved)
            time_array = time_array + (delta_sec - t0_delay).value * SECPERDAY

            # replace time array <- wasn't being restuffed before
            in_hdul['EVENTS', e_indx].data[in_col] = time_array

            # add delaytime to EXPSTART and EXPEND, and update header
            for keyword in ['EXPSTART', 'EXPEND']:
                mjd = events_tab.header[keyword]
                if hst_orb is None:
                    delta_sec = calc_delay_jpl(mjd, ra, dec, distance=distance)
                else:
                    delta_sec = calc_delay_orbfile(mjd, ra, dec, hst_orb, distance=distance)

                # If converting to TDB, apply clock correction. Otherwise, keep as is.
                if time_system == 'UTC':
                    events_tab.header[keyword] = mjd + delta_sec.value
                elif time_system == 'TDB':
                    mjd_obj = Time(mjd, format='mjd', scale='utc')
                    events_tab.header[keyword] = mjd_obj.tdb.value + delta_sec.value

            in_hdul.flush()
            if verbose:
                print(f"    [EVENTS,{e_indx}] extension has been updated")

        # Loop through all x1d or image extensions (there will be none if
        # the input is an events file).  Note that we only expect EXPSTART
        # and EXPEND to be present in SCI extensions; however, they
        # could be in ERR and DQ as well (e.g. in output from inttag),
        # and if they are present they must be updated.
        for e_indx in range(1, nsci_ext + 1):
            cur_tab = in_hdul[e_indx]

            modified = False

            # add delaytime to EXPSTART and EXPEND, and update header
            for keyword in ['EXPSTART', 'EXPEND']:
                if keyword in cur_tab.header:
                    mjd = cur_tab.header[keyword]
                    if hst_orb is None:
                        delta_sec = calc_delay_jpl(mjd, ra, dec, distance=distance)
                    else:
                        delta_sec = calc_delay_orbfile(mjd, ra, dec, hst_orb, distance=distance)

                    # If converting to TDB, apply clock correction. Otherwise, keep as is.
                    if time_system == 'UTC':
                        cur_tab.header[keyword] = mjd + delta_sec.value
                    elif time_system == 'TDB':
                        mjd_obj = Time(mjd, format='mjd', scale='utc')
                        cur_tab.header[keyword] = mjd_obj.tdb.value + delta_sec.value

                    modified = True

            in_hdul.flush()
            if verbose and modified:
                print(f"    extension {e_indx} has been updated")

        if time_script:
            tcheck2 = time.time()
            print(f'Checkpoint 2: {tcheck2 - tstart} s')

        # IF STIS
        # add delaytime to TEXPSTRT and TEXPEND, and update primary header
        # COS has no TEXPSTRT in primary header
        if in_hdul[0].header['INSTRUME'] == "STIS":
            for keyword in ['TEXPSTRT', 'TEXPEND']:
                mjd = in_hdul[0].header[keyword]
                if hst_orb is None:
                    delta_sec = calc_delay_jpl(mjd, ra, dec, distance=distance)
                else:
                    delta_sec = calc_delay_orbfile(mjd, ra, dec, hst_orb, distance=distance)

                # If converting to TDB, apply clock correction. Otherwise, keep as is.
                if time_system == 'UTC':
                    in_hdul[0].header[keyword] = mjd + delta_sec.value
                elif time_system == 'TDB':
                    mjd_obj = Time(mjd, format='mjd', scale='utc')
                    in_hdul[0].header[keyword] = mjd_obj.tdb.value + delta_sec.value

        # add keyword to flag the fact that the times have been corrected
        in_hdul[0].header.insert('DQICORR', ('DLAYCORR', 'COMPLETE', 'delaytime has been applied'), after=True)
        in_hdul[0].header['history'] = "Times corrected to solar system barycenter;"
        in_hdul[0].header['history'] = f"All times now in BJD_{time_system};"

        if hst_orb is not None:
            in_hdul[0].header['history'] = f"HST orbital table {hst_orb}"
        else:
            in_hdul[0].header['history'] = "no HST orbital tables were used."

        # close in file
        in_hdul.close()

        if verbose:
            print("... done")

        if time_script:
            tcheck3 = time.time()
            print(f'Checkpoint 3: {tcheck3 - tstart} s')


def calc_delay_jpl(times, ra, dec, distance=1e9, verbose=True):
    r"""
    Calculate the barycentric light-travel time correction for HST using JPL Horizons.

    Compute the barycentric light-travel time delay for the Hubble Space Telescope (HST)
    at the provided observation times and a target sky position (RA, Dec) by querying
    JPL Horizons for HST's geocentric state-vector (or interpolating a regular-sampled
    vector set) and combining it with Earth's barycentric position. A finite-distance
    correction is applied to account for targets that are not at infinite distance.

    Parameters
    ----------
    times: array-like or `~astropy.time.Time`
        Observation times in Modified Julian Date (MJD). Can be a scalar or an array.

    ra: float or `~astropy.units.Quantity`
        Right ascension of the target in degrees.

    dec: float or `~astropy.units.Quantity`
        Declination of the target in degrees.

    distance: float, optional
        Distance to the target (default ``1e9``). This value is used in the finite-distance
        correction term. (See ``Notes`` for how it enters the equation.)

    verbose: bool, optional
        If True (default) print diagnostics about interpolation, the finite-distance
        correction, and the computed light-travel times.

    Returns
    -------
    lt_time: `~astropy.units.Quantity`
        Light-travel time correction(s) (Astropy Quantity) in days. This is the time
        that should be added to the observation times to obtain barycentric arrival times,
        and includes the finite-distance correction.

    Notes
    -----
    - Uses the JPL planetary ephemeris (``jpl``) via Astropy's ``solar_system_ephemeris``
      for consistency with Horizons.
    - For multiple times the function queries Horizons for HST (-48, location 500)
      with a regular step (1 minute) over the time range (plus a 5-minute margin)
      and cubic-interpolates the returned vectors to the requested times. For a single
      scalar time it queries Horizons directly at that epoch.
    - The finite-distance correction implemented here follows the same algebraic
      form used elsewhere in this codebase:

      .. math::

           \mbox{correction\_term} = \frac{-0.5}{D} \left(|r|^2 - (n \cdot r)^2\right) / c

      where `r` is the HST barycentric position vector, `n` is the unit vector toward
      the target, `D` is the provided `distance`, and `c` is the speed of light.
      The expression is converted into days before being returned.
    - The HST vector returned from Horizons is treated as geocentric and then transformed
      to ITRS/ICRS and added to Earth's barycentric position from Astropy's
      ``get_body_barycentric('earth', ...)`` to obtain HST's barycentric position.
    - Requires `astroquery` (Horizons), Astropy with the `jplephem` ephemeris available,
      and a working internet connection for Horizons queries when not using a local orbit file.

    **Performance / Accuracy Considerations**

    - Interpolation at 1-minute sampling is a compromise between speed and accuracy;
      a timing error of ~1 minute corresponds to a maximum light-travel time error of
      a few milliseconds. Increase sampling density or query exact epochs if sub-ms
      accuracy is required.
    - Converting large arrays of epochs and coordinate transforms may be time-consuming.

    Raises
    ------
    ValueError
        If the Horizons query or interpolation fails to produce vectors covering the
        requested time range, or if the vector transformation cannot be performed.

    RemoteServiceError
        If the Horizons query or interpolation fails to produce vectors covering the
        requested time range.

    ConvertError
        If the vector transformation cannot be performed.

    Examples
    --------
    >>> # single time (MJD), RA/Dec in degrees
    >>> barycentric_correction.calc_delay_jpl(55521.123, 210.8023, -47.393)
    Finite distance correction: -7.804152203899432e-08 s
    Light travel times: -405.56807173759836 s

    >>> # multiple times
    >>> times = [60200.123, 60200.124, 60200.125]
    >>> lt_array = calc_delay_jpl(times, 210.8023, -47.393, distance=1e9, verbose=False)
    """

    # Using the JPL epehermis to be consistent with what Horizons gives
    # see https://github.com/astropy/astropy/pull/11608
    # Will require jplephem package, but that's already in stenv!
    solar_system_ephemeris.set('jpl')

    # Put times into Astropy time object
    # We will query HST's position at this time
    # We will be off by the HST-geocenter light travel time, FWIW
    times_geo = Time(times,
                     format='mjd',
                     scale='utc',
                     location=EarthLocation.from_geocentric(0, 0, 0, unit='m'))

    # Get HST's location

    # Interpolate HST's position if more than just one time:
    # We can't query all times to Horizons because there's
    # too many. Let's get HST's position every minutes instead,
    # and then interpolate...
    # Position every minute should be good enough. If we were
    # off by a minute, that would be a max of 6ms light travel time
    # So with interpolation, we'll be dandy.
    # Adding one minute to each side to avoid extrapolating

    # Interpolate HST's position if more than just one time
    if np.isscalar(times) is False:
        if verbose:
            print('Multiple times: interpolating positions')

        # Query five minute before and after min and max times
        # (In case we only query a minute of tags, we need more for spline)
        # With a step of 1 min
        epochs = {'start': (np.min(times_geo.tdb) - (5.0 * u.min)).
                  to_datetime().strftime("%Y-%m-%d %H:%M:%S"),
                  'stop': (np.max(times_geo.tdb) + (5.0 * u.min)).
                  to_datetime().strftime("%Y-%m-%d %H:%M:%S"),
                  'step': '1m'}
        hstobj = Horizons(id='-48', location='500', epochs=epochs)

        # Must be refplane = 'earth' not 'ecliptic
        hstvec = hstobj.vectors(refplane='earth')

        # Run interpolation
        f_hstx = interp1d(hstvec['datetime_jd'], hstvec['x'], kind='cubic')
        f_hsty = interp1d(hstvec['datetime_jd'], hstvec['y'], kind='cubic')
        f_hstz = interp1d(hstvec['datetime_jd'], hstvec['z'], kind='cubic')

        hstvecx = f_hstx(times_geo.tdb.jd)
        hstvecy = f_hsty(times_geo.tdb.jd)
        hstvecz = f_hstz(times_geo.tdb.jd)

        # Re-package with units
        hstarr = [hstvecx, hstvecy, hstvecz] * u.AU

    else:
        # Query HST's gencentric position in JPL Horizons
        hstobj = Horizons(id='-48', location='500', epochs=times_geo.tdb.jd)

        # Must be refplane = 'earth' not 'ecliptic
        hstvec = hstobj.vectors(refplane='earth')

        # Re-package with units
        hstarr = [hstvec['x'][0], hstvec['y'][0], hstvec['z'][0]] * u.AU

    # Convert from vector table to Astropy Quantity
    # Then convert from GCRS to ITRS
    # i.e., from Earth-centric to Barycentric
    # This will still be the distance from HST to Earth,
    # but in ICRS directions
    # hstarr = [hstvec['x'][0], hstvec['y'][0], hstvec['z'][0]]*u.AU

    with erfa_astrom.set(ErfaAstromInterpolator(1000 * u.s)):

        # Make into GCRS object at proper times
        hstGCRS = GCRS([hstarr[0],
                        hstarr[1],
                        hstarr[2]],
                       representation_type='cartesian',
                       obstime=times_geo)

        # Convert to ICRS so we can add to Earth's barycentric location
        # This takes a bit of time if there are many times to convert
        hstITRS = hstGCRS.transform_to(ITRS(obstime=times_geo))

        # Get Earth's position
        # Can also take a few minutes with lots of times
        earthICRS = get_body_barycentric('earth', times_geo)

        # Add the two vectors together to get HST's position relative
        # to the Solar System barycenter
        hstbary = [earthICRS.x.to('AU').value + hstITRS.x.value,
                   earthICRS.y.to('AU').value + hstITRS.y.value,
                   earthICRS.z.to('AU').value + hstITRS.z.value]

        # Define our targets location on the sky
        target = SkyCoord(ra, dec,
                          unit=(u.deg, u.deg), frame='icrs')

        # Put into an array to dot with hstbary
        # Cartesian makes this a unit vector in direction of target
        target_arr = [target.cartesian.x.value,
                      target.cartesian.y.value,
                      target.cartesian.z.value]

        # Calculate the finite-distance correction term
        # make the sum over the correction axis
        correction_term = ((-0.5 / distance) *
                           (np.sum((np.array(hstbary))**2, axis=0) -
                                  (np.dot(target_arr, hstbary))**2) * u.AU /
                           c.c).to('day')

        if verbose:
            if correction_term.size < 10:
                print(f"Finite distance correction: {correction_term.to('s')}")
            else:
                print(f"Finite distance correction: \
                      {correction_term[0:10].to('s')}")

        # Let's now define HST's location relative to the barycenter
        hstloc = EarthLocation.from_geocentric(x=hstITRS.x,
                                               y=hstITRS.y,
                                               z=hstITRS.z)
        # Define the times, now with the correct location
        hsttime = Time(times, format='mjd', scale='utc', location=hstloc)

        # Calculate the light travel time,
        # adding the correction term above.
        # Then define the new barycenter times!
        lt_time = hsttime.light_travel_time(target) + correction_term

    if verbose:
        if lt_time.size < 10:
            print(f"Light travel times: {lt_time.to('s')}")
        else:
            print(f"Light travel times: {lt_time[0:10].to('s')}")

    return lt_time


def calc_delay_orbfile(times, ra, dec, hst_orb, distance=1e9, verbose=True):
    r"""
    Calculate the light-travel time correction for HST observations using an orbit file.

    This function computes the barycentric light-travel time delay for the Hubble Space Telescope (HST)
    given a set of observation times and a target sky position (RA, Dec). It interpolates HST's position
    from a provided orbit file, combines it with Earth's barycentric position, and calculates the
    light-travel time correction to the Solar System barycenter, including a finite-distance correction term.

    Parameters
    ----------
    times: array-like or float
        Observation times in Modified Julian Date (MJD), corresponding to HST exposures.

    ra: float
        Right ascension of the target in degrees.

    dec: float
        Declination of the target in degrees.

    hst_orb: str
        Path to the HST orbit FITS file. This file must contain columns `TIME`, `X`, `Y`, and `Z`
        giving HST's position (in km) relative to the Earth's center.

    distance: float, optional
        Distance to the target in kilometers (default is `1e9`, effectively infinite distance).
        Used to apply the finite-distance light-travel time correction.

    verbose: bool, optional
        If True (default), print information about the finite-distance correction
        and the calculated light-travel times.

    Returns
    -------
    lt_time: `~astropy.units.Quantity`
        The barycentric light-travel time correction(s) in days, including the finite-distance correction term.

    Notes
    -----
    - Uses the JPL planetary ephemeris (`jplephem`) for consistency with NASA Horizons.
    - Interpolates HST's orbital position at each observation time using cubic interpolation
      to avoid excessive Horizons queries.
    - The resulting correction term accounts for HST's motion relative to the Solar System barycenter.
    - The finite-distance correction term scales as ``~(r^2 - (r·n)^2) / (2cD)``, where
      `r` is HST's barycentric position vector, `n` is the unit vector toward the target,
      `c` is the speed of light, and `D` is the target distance.
    - HST orbit files may be downloaded using information provided at
      https://www.stsci.edu/~STIS/monitors/ephemeris_files.html

    Raises
    ------
    OrbFileError
        If the orbit file does not cover the requested time range.

    Examples
    --------
    >>> barycentric_correction.calc_delay_orbfile([55521.123], 210.8023, -47.393, hst_orb='pubj0000r.fit')
    Finite distance correction: [-7.80563708e-08] s
    Light travel times: [-405.56687377] s
    """

    if not os.access(hst_orb, os.F_OK):
        raise FileNotFoundError('HST orbfile not found')

    # Using the JPL epehermis to be consistent with what Horizons gives
    # see https://github.com/astropy/astropy/pull/11608
    # Will require jplephem package, but that's already in stenv!
    solar_system_ephemeris.set('jpl')

    # Put times into Astropy time object
    # We will query HST's position at this time
    # We will be off by the HST-geocenter light travel time, FWIW
    times_geo = Time(times,
                     format='mjd',
                     scale='utc',
                     location=EarthLocation.from_geocentric(0, 0, 0, unit='m'))

    # Get HST's location

    # Interpolate HST's position if more than just one time:
    # We can't query all times to Horizons because there's
    # too many. Let's get HST's position every minutes instead,
    # and then interpolate...
    # Position every minute should be good enough. If we were
    # off by a minute, that would be a max of 6ms light travel time
    # So with interpolation, we'll be dandy.
    # Adding one minute to each side to avoid extrapolating

    # HST orb file
    with fits.open(hst_orb) as hdu_orb:
        # Handle case changes in orbital file time column name:
        try:
            in_col = [x.name for x in hdu_orb[1].data.columns if x.name.lower().strip() == 'time'][0]
        except IndexError as e:
            raise ValueError(f'Unable to locate "Time" column in hst_orb file "{hst_orb}".') from e

        f_hstx = interp1d(float(hdu_orb[1].header['FIRSTMJD']) + hdu_orb[1].data[in_col].astype(np.float64) / SECPERDAY, hdu_orb[1].data['X'], kind='cubic')
        f_hsty = interp1d(float(hdu_orb[1].header['FIRSTMJD']) + hdu_orb[1].data[in_col].astype(np.float64) / SECPERDAY, hdu_orb[1].data['Y'], kind='cubic')
        f_hstz = interp1d(float(hdu_orb[1].header['FIRSTMJD']) + hdu_orb[1].data[in_col].astype(np.float64) / SECPERDAY, hdu_orb[1].data['Z'], kind='cubic')

        try:
            hstvecx = f_hstx(times_geo.tdb.mjd)
            hstvecy = f_hsty(times_geo.tdb.mjd)
            hstvecz = f_hstz(times_geo.tdb.mjd)
        except ValueError as e:
            raise OrbFileError('Orbit file does not match observation times. '
                               'Make sure you have got the correct one.') from e

        hstarr = [hstvecx, hstvecy, hstvecz] * u.km

    # Convert from vector table to Astropy Quantity
    # Then convert from GCRS to ITRS
    # i.e., from Earth-centric to Barycentric
    # This will still be the distance from HST to Earth,
    # but in ITRS directions
    # hstarr = [hstvec['x'][0], hstvec['y'][0], hstvec['z'][0]] * u.AU

    with erfa_astrom.set(ErfaAstromInterpolator(1000 * u.s)):

        # Make into ICRS object at proper times
        hstICRS = ICRS([hstarr[0],
                        hstarr[1],
                        hstarr[2]],
                       representation_type='cartesian')

        # Get Earth's position
        # Can also take a few minutes with lots of times
        earthICRS = get_body_barycentric('earth', times_geo)

        # Add the two vectors together to get HST's position relative
        # to the Solar System barycenter
        # We keep this in an array format to calculate the correction term,
        # but convert to ICRS coordinate for light travel times.
        hstbary = [earthICRS.x.to('AU').value + hstICRS.x.to('AU').value,
                   earthICRS.y.to('AU').value + hstICRS.y.to('AU').value,
                   earthICRS.z.to('AU').value + hstICRS.z.to('AU').value]

        hstbary_icrs = ICRS(hstbary * u.AU, representation_type='cartesian')

        # from_geocentric() is expecting an ITRS coordinate, so we need to convert.
        hstbary_itrs = hstbary_icrs.transform_to(ITRS(obstime=times_geo))
        itrs_cartesian = hstbary_itrs.cartesian

        g = EarthLocation.from_geocentric(x=itrs_cartesian.x,
                                          y=itrs_cartesian.y,
                                          z=itrs_cartesian.z)

        # Define our targets location on the sky
        target = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')

        # Put into an array to dot with hstbary
        # Cartesian makes this a unit vector in direction of target
        target_arr = [target.cartesian.x.value,
                      target.cartesian.y.value,
                      target.cartesian.z.value]

        # Calculate the finite-distance correction term
        correction_term = ((-0.5 / distance) *
                           (np.sum((np.array(hstbary))**2) -
                            (np.dot(target_arr, hstbary))**2) * u.AU /
                           c.c).to('day')

        if verbose:
            if correction_term.size < 10:
                print(f"Finite distance correction: {correction_term.to('s')}")
            else:
                print(f"Finite distance correction: \
                      {correction_term[:10].to('s')}")

        # Let's now define HST's location relative to the barycenter
        # from_geocentric() is expecting an ITRS coordinate!
        hstloc = EarthLocation.from_geocentric(x=itrs_cartesian.x,
                                               y=itrs_cartesian.y,
                                               z=itrs_cartesian.z)
        # Define the times, now with the correct location
        hsttime = Time(times, format='mjd', scale='utc', location=hstloc)

        # Calculate the light travel time,
        # adding the correction term above.
        # Then define the new barycenter times!
        lt_time = hsttime.light_travel_time(target) + correction_term

    if verbose:
        if lt_time.size < 10:
            print(f"Light travel times: {lt_time.to('s')}")
        else:
            print(f"Light travel times: {lt_time[0:10].to('s')}")

    return lt_time


def combine_hst_orb_files(file1, file2, outname, overwrite=False):
    """
    Combine two neighboring hst_orb ephemeris files (ORBs) together.  Useful for
    correcting MAST datasets that are observed over an ORB file boundary.

    Parameters
    ----------
    file1: str
        Path to the first hst_orb FITS file.

    file2: str
        Path to the second hst_orb FITS file.

    outname: str
        Name/path of combined hst_orb FITS file to output.

    overwrite: bool
        Overwrite a pre-existing output file?  Default=False

    Returns
    -------
    None

    Raises
    ------
    OrbFileError
        The time gap between orbit files is too large (>120 s)

    FileExistsError
        Output file already exists and overwrite=False

    ValueError
        Unexpected duplicate times detected after truncation of later file
    """

    if not overwrite and os.access(outname, os.F_OK):
        raise FileExistsError(f"Output file already exists.  Remove or set overwrite=True.")

    # Stack arrays together, excluding the repeated point in the middle:
    t1 = Table.read(file1, hdu=1)
    t2 = Table.read(file2, hdu=1)

    # Check time ordering and contiguousness of files:
    if t1.meta['FIRSTMJD'] > t2.meta['FIRSTMJD']:
        tmpfile = file1
        file1 = file2
        file2 = tmpfile

        tmp = t1
        t1 = t2
        t2 = tmp

        del tmpfile, tmp

    acceptable_gap = 120.  # s
    if (t2.meta['FIRSTMJD'] - t1.meta['LASTMJD']) > acceptable_gap / SECPERDAY:
        raise OrbFileError(f"hst_orb files have a gap > {acceptable_gap:.0f} s")

    # Access time arrays case-insensitively:
    time_col = [x for x in t1.columns if x.lower().strip() == 'time'][0]
    time_col2 = [x for x in t2.columns if x.lower().strip() == 'time'][0]
    t2.rename_column(time_col2, time_col)

    # Increment later file's time array to be relative to the earlier file:
    t2[time_col] += (t2.meta['FIRSTMJD'] - t1.meta['FIRSTMJD']) * SECPERDAY

    # Concatenate non-overlapping times:
    g = t2[time_col] > t1[time_col].max()
    # Silence MergeConflictWarnings of unused output metadata:
    with warnings.catch_warnings(record=True) as _:
        t_combined = vstack([t1, t2[g]])
    if len(t_combined) != len(set(t_combined[time_col])):
        raise ValueError('Non-unique times found in input hst_orb files')

    # Combine into a single output FITS file:
    hdu = fits.HDUList()
    with fits.open(file1) as f1, fits.open(file2) as f2:
        hdr0 = f1[0].header.copy()
        hdr1 = f1[1].header.copy()

        # Update metadata:
        hdr0['FILENAME'] = os.path.basename(outname)
        hdr0['DATE'] = datetime.datetime.now().isoformat().rsplit('.', 1)[0]
        hdr0['LEAPSECO'] = f1[0].header.get('LEAPSECO', False) and f2[0].header.get('LEAPSECO', False)
        hdr0['LEAPSECS'] = f1[0].header.get('LEAPSECS', 0) + f2[0].header.get('LEAPSECS', 0)
        hdr1['ROOTNAME'] = f"{f1[1].header.get('ROOTNAME', 'UNKNOWN').strip()} + {f2[1].header.get('ROOTNAME', 'UNKNOWN').strip()}"
        hdr1['EXPNAME'] = f"{f1[1].header.get('EXPNAME', 'UNKNOWN').strip()} + {f2[1].header.get('EXPNAME', 'UNKNOWN').strip()}"
        # Update keywords in earlier file's header with appropriate end times:
        hdr0['DEFINEND'] = f2[0].header.get('DEFINEND', 'UNKNOWN')
        hdr1['LASTMJD'] = f2[1].header.get('LASTMJD', 'UNKNOWN')

        hdr0.add_history('This file was created from components via')
        hdr0.add_history('stistools.barycentric_correction.combine_hst_orb_files() from:')
        hdr0.add_history(os.path.basename(file1))
        hdr0.add_history(os.path.basename(file2))

        hdu.append(fits.PrimaryHDU(header=hdr0))
        hdu.append(fits.BinTableHDU(data=t_combined, header=hdr1))
        hdu.writeto(outname, checksum=True, overwrite=True)

    print(f"Combined hst_orb file written:  {outname}")


def odelay_file_compare(file1, file2):
    """
    Compare timing information between two FITS files.

    Computes and prints the differences in exposure start times and data timestamps
    between two FITS files, typically used for verifying time coordinate consistency
    in astronomical observations.

    Parameters
    ----------
    file1: str
        Path to the first FITS file.

    file2: str
        Path to the second FITS file to compare against file1.

    Returns
    -------
    None
        This function only prints comparison results to stdout.

    Notes
    -----
    The function prints the following comparisons (file2 - file1):

    - TEXPSTRT header keyword difference in days and seconds
    - EXPSTART header keyword difference in days and seconds (from extension 1)
    - First and last TIME values difference in seconds (only if 'tag' is in file1 name)
    - TT to TDB time scale conversion difference for reference

    The function assumes:

    - Both files have TEXPSTRT in the primary header (extension 0)
    - Both files have EXPSTART in extension 1 header
    - If ``tag`` appears in file1 name, extension 1 contains a table.
      The script will print the difference between the first and
      last times in each file.

    Examples
    --------
    ::

        >>> odelay_file_compare('observation1_tag.fits', 'observation2_tag.fits')
        TEXPSTRT f2-f1: -1.5854311641305685e-08 days
        TEXPSTRT f2-f1: -0.0013698125258088112 seconds
        EXPSTART f2-f1: -1.5854311641305685e-08 days
        EXPSTART f2-f1: -0.0013698125258088112 seconds
        First TIME f2-f1: 0.0 seconds
        Last TIME f2-f1: -0.23225000000002183 seconds
        In case it is helpful, the difference between TT and TDB_BJD is -0.001497 s
    """

    with fits.open(file1) as f1, fits.open(file2) as f2:
        # IF STIS
        # add delaytime to TEXPSTRT and TEXPEND, and update primary header
        # COS has no TEXPSTRT in primary header
        if f1[0].header['INSTRUME'] == "STIS":
            print(f"TEXPSTRT file2-file1 {(f2[0].header['TEXPSTRT'] - f1[0].header['TEXPSTRT']) * SECPERDAY} seconds")

        print(f"EXPSTART file2-file1: {(f2[1].header['EXPSTART'] - f1[1].header['EXPSTART']) * SECPERDAY} seconds")

        if 'tag' in os.path.basename(file1):
            # Handle case changes in file's time column name:
            try:
                in_col = [x.name for x in f1[1].data.columns if x.name.lower().strip() == 'time'][0]
            except IndexError as e:
                raise ValueError(f'Unable to locate "TIME" column in file1 "{file1}".') from e

            print(f"First TIME file2-file1: {f2[1].data[in_col][0] - f1[1].data[in_col][0]} seconds")
            print(f"Last TIME file2-file1: {f2[1].data[in_col][-1] - f1[1].data[in_col][-1]} seconds")

        if f1[0].header['INSTRUME'] == "STIS":
            t = Time(f1[0].header['TEXPSTRT'], format='mjd', scale='utc')
        elif f1[0].header['INSTRUME'] == "COS":
            t = Time(f1[1].header['EXPSTART'], format='mjd', scale='utc')
        else:
            raise ValueError(f"Unexpected INSTRUME value: {f1[1].header['EXPSTART']}")

    diff = (t.tdb.value - t.tt.value) * SECPERDAY
    print(f'In case it is helpful, the difference between TT and TDB_BJD is {diff:.6f} s')
