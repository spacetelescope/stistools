#! /usr/bin/env python

from datetime import datetime
from math import cos, sin
import numpy as np
from scipy.interpolate import lagrange, interp1d

from astropy.io import fits
from astropy.table import Table, vstack, hstack
from astropy.units import cds

# Lots of notes in the original file

'''
Compute the times of observation relative to the solar-system barycenter
of a distant object (outside the solar system).  It is assumed that
the set of observations covers a short enough time that the object does
not move appreciably during the observations.  Note that a one-arcsecond
error in the position of the object can result in a time-delay error
of ~2.4 millisec.  It is also assumed that the set of observations covers
a short enough time interval that the relativistic correction, if required,
remains constant during the observations.

#  Date		Author		Description
#  ----		------		-----------
#  10-Nov-1984	C. D. Biemesderfer	Original module
#  18-Apr-1990  J.-C. Hsu	rewrite in SPP
#  19-Jun-1992  J.-C. Hsu	read RA_TARG and DEC_TARG from header
#  19-Aug-1997  J.-C. Hsu	modify for STIS data
#  23-Aug-2000  Phil Hodge	remove out_col & add verbose arguments;
#				modify times in-place; update GTI table and
#				header keywords as well; change INTERVAL
#				from 1. to 0.1; add DELAYCOR and HISTORY
#  10-Jun-2003  Phil Hodge	modify to work with SCI extensions;
#				change calling sequence (ephem table names
#				instead of table pointers)



pointer	fin			# file template pointer
int	nfin 			# number of files in the input template
double	parallax		# parallax of the target (in arc sec)
char	earth_ephem[SZ_FNAME]	# earth ephemeris table name
char	obs_ephem[SZ_FNAME]	# list of observer ephemeris table names
char	in_col[SZ_LINE]		# time column name
bool	verbose			# print timing info?


pointer	fin			# i: file template pointer
int	nfin 			# i: number of files in the input template
double	parallax		# i: parallax of the target (in arc sec)
char	earth_ephem[ARB]	# i: name of the earth_ephem table
char	obs_ephem[ARB]		# i: names of the obs_ephem tables
char	in_col[ARB]		# i: time column name
#--
pointer sp
pointer history		# for constructing a history record
pointer extname		# from first extension, for checking file type
double	ra, dec		# RA and Dec of the target at J2000 (in deg.)
double	objvec[3]	# unit geocentric state vector of the target
double	mjd1, mjd2 	# MJD of the observation start/end time
double	epoch		# MJD of an event
double	delta_sec	# correction (sec) to be added to TIME col
double	t0_delay	# correction at TEXPSTRT
char	ifile[SZ_FNAME]
char	ifilen[SZ_FNAME]
int	nchar, i, j, k
pointer	time_e, x_e, y_e, z_e, time_o, x_o, y_o, z_o
pointer	tp		# pointer of the input table
pointer	cptr		# for the TIME column in EVENTS tables
pointer cp1, cp2	# for START and STOP columns in GTI table
int	nrows		# number of rows in a table
int	nextend		# keyword from primary header
int	nevents_tab	# number of EVENTS tables (one less than NEXTEND, or 0)
int	nsci_ext	# number of x1d or image extensions (0 or NEXTEND)
int	npts_earth, npts_obs
double	tm		# a time from the TIME column
double	tm_prev		# the last tm for which all_delay was called
char	text[SZ_TIME]	# for printing timing info
char	delaycorr[SZ_FNAME]
int	frac		# for printing percentage done
int	filetype	# events table, x1d table, image set
bool	modified	# true if an x1d or image extension was updated



'''

NLAN_EARTH = 10
NLAN_OBS = 10
INTERVAL = 0.1      # in seconds

# File types, probably don't need this...
EVENTS_TABLE = 1

SECPERDAY = 86400  #D0	# number of sec in a day
MINPERDAY = 1440   #D0		# number of min in a day
HRPERDAY = 24     #D0		# number of hours in a day

LYPERPC	= 3.261633   #D0	# light years per parsec
KMPERPC	= 3.085678e13  #D13	# kilometers per parsec
AUPERPC	= 206264.8062470964  #D0	# how many AU in one parsec
KMPERAU	= KMPERPC/AUPERPC
CLIGHT = 499.00479   #D0	# light travel time (in sec) of 1 AU

# might not need these
EARTH_EPHEMERIS = 1
OBS_EPHEMERIS = 2

JD_TO_MJD = 2400000.5   #d0	# subtract from JD to get MJD

RADIAN = 57.295779513082320877

# Need to figure out what/if default for distance parameter

def odelaytime(table_names, earth_ephem, obs_ephem=None, distance=100000.0,
               dist_unit="pc", in_col="TIME", verbose=False):

    parallax = distance

    # convert distance to parsec first
    if dist_unit == "arcsec":
        pass
    elif dist_unit == "au":
        parallax /= AUPERPC
    elif dist_unit == "ly":
        parallax /= LYPERPC
    elif dist_unit == "km":
        parallax /= KMPERPC
    else:
        if dist_unit != "pc":
            raise ValueError("illegal distance unit: {}".format(dist_unit))

    # check that distance is non-zero (pc)
    if parallax <= 0:
        raise ValueError("non-positive distance: {}".format(parallax))

    # I don't think I need this part
    # call smark(sp)
    # call salloc(history, SZ_FNAME, TY_CHAR)
    # call salloc(extname, SZ_FNAME, TY_CHAR)

    # time_e, x_e, y_e, z_e, npts_earth)
    earth_ephem_table = get_ephem(earth_ephem, EARTH_EPHEMERIS)
    # time_o, x_o, y_o, z_o, npts_obs
    # obs_ephem_table will be None if no file names were provided
    obs_ephem_table = get_ephem(obs_ephem, OBS_EPHEMERIS)

    if obs_ephem_table is None:
        npts_obs = 0
    else:
        npts_obs = len(obs_ephem_table)

    for in_table_file in table_names:
        if verbose:
            print("odelaytime: processing {} ...".format(in_table_file))

        # determine the file type, based on the first extension
        # original code (ln 124) contains some kind of catch all
        # for SCI_IMAGE? It also caught anything outside of those options,
        # might want to add that back in as well.
        in_hdul = fits.open(in_table_file, mode='update')
        extname = in_hdul[1].header['EXTNAME']
        if extname == "EVENTS":
            filetype = "EVENTS_TABLE"
        elif extname == "SCI":
            filetype = "X1D_TABLE"

        if 'DELAYCOR' in in_hdul[0].header and \
                        in_hdul[0].header['DELAYCOR'] == "COMPLETE":
            print("{} has already been corrected, so no further processing"
                  " will be applied to this file".format(in_table_file))
            continue
        else:
            in_hdul[0].header['DELAYCOR'] = "PERFORM"

        mjd1 = in_hdul[0].header['TEXPSTRT']
        mjd2 = in_hdul[0].header['TEXPEND']
        ra = in_hdul[0].header['RA_TARG']
        dec = in_hdul[0].header['DEC_TARG']

        if 'NEXTEND' in in_hdul[0].header:
            nextend = in_hdul[0].header['NEXTEND']
        else:
            nextend = 1

        if filetype == "EVENTS_TABLE":
            # assume (for now) that the last extension is a GTI table
            nevents_tab = max(nextend-1, 1)
            nsci_ext = 0
        elif filetype == "X1D_TYPE":
            nevents_tab = 0
            nsci_ext = nextend
        else:
            nevents_tab = 0
            nsci_ext = nextend

        # check the range
        if mjd1 < earth_ephem_table['Time'][0] or \
                        mjd2 > earth_ephem_table['Time'][-1]:
            print("Epoch is outside the Earth ephemeris range so {} will be"
                  " skipped".format(in_table_file))
            continue

        # this checks npts_obs first
        if obs_ephem_table is not None:
            if mjd1 < obs_ephem_table['Time'][0] or \
                            mjd2 > obs_ephem_table['Time'][-1]:
                print("Epoch is outside the obs ephemeris range so {} will be"
                      " skipped".format(in_table_file))


        # Check that all expected extensions and the time column are
        # present.  We do this here before making any change, to avoid
        # the possibility of correcting only part of the input file.
        if not ext_exist(in_table_file, in_col, nevents_tab, nsci_ext):
            print("missing extensions... ")
            continue

        # calculate the target's positional vector
        cosdec = cos(dec / RADIAN)
        objvec = np.array([cosdec * cos(ra / RADIAN),
                            cosdec * sin(ra / RADIAN),
                            sin(dec / RADIAN)])

        # Get the delay at the exposure start time.  We'll subtract this
        # delay from each time that we update, if that time is relative
        # to EXPSTART (which must be the same as TEXPSTRT).
        t0_delay = all_delay(mjd1, parallax, objvec, earth_ephem_table,
                             len(earth_ephem_table), obs_ephem_table, npts_obs)

        if filetype == "EVENTS_TABLE":
            # update times in GTI table
            if 'GTI' in in_hdul:
                gti_tab = in_hdul['GTI'].data
                nrows = len(gti_tab)
                for row in gti_tab:
                    tm = row['START']
                    epoch = mjd1 + tm / SECPERDAY
                    delta_sec = all_delay(epoch, parallax, objvec,
                                          earth_ephem_table,
                                          len(earth_ephem_table),
                                          obs_ephem_table,
                                          npts_obs)
                    row['START'] = tm + (delta_sec - t0_delay)

                    tm = row['STOP']
                    epoch = mjd1 + tm / SECPERDAY
                    delta_sec = all_delay(epoch, parallax, objvec,
                                          earth_ephem_table,
                                          len(earth_ephem_table),
                                          obs_ephem_table,
                                          npts_obs)
                    row['STOP'] = tm + (delta_sec - t0_delay)

                in_hdul.flush()
                if verbose:
                    print("  GIT extension has been updated")

            else:
                # If there's no GTI, last extensions is an EVENTS table
                nevents_tab += 1

        # loop through all EVENTS extensions (there will be none if
        # the input is an x1d or image file)
        for e_indx in range(1, nevents_tab+1):
            events_tab = in_hdul['EVENTS', e_indx]
            mjd1 = events_tab.header['EXPSTART']
            mjd2 = events_tab.header['EXPEND']

            # find the time (since EXPSTART) column in the table
            nrows = len(events_tab.data[in_col])

            tm_prev = 0

            # print out timing info
            if verbose:
                print("  start time: ", datetime.now())

            # pull out just the time array
            time_array = events_tab.data['TIME']

            # go through each row
            print(nrows)
            for r_indx in range(nrows):
                tm = time_array[r_indx]

                # if the time is within INTERVAL of the previous time,
                # apply the same delaytime, otherwise recalculate
                if r_indx == 0 or (tm - tm_prev) > INTERVAL:
                    epoch = mjd1 + tm / SECPERDAY

                    delta_sec = all_delay(epoch, parallax, objvec,
                                          earth_ephem_table,
                                          len(earth_ephem_table),
                                          obs_ephem_table,
                                          npts_obs)
                    tm_prev = tm

                # write the corrected time back to the time column
                time_array[r_indx] = tm + (delta_sec - t0_delay)

                if verbose and r_indx in np.arange(0, 100000, 10000):
                    # print percent done
                    # print("    Percentage done:  {}".
                    #      format(int(100*(r_indx/nrows))))
                    #print(r_indx)
                    pass

            # re-stuff time array
            in_hdul['EVENTS', e_indx].data[in_col] = time_array

            # add delaytime to EXPSTART and EXPEND, and update header
            delta_sec = all_delay(mjd1, parallax, objvec,
                                  earth_ephem_table,
                                  len(earth_ephem_table),
                                  obs_ephem_table,
                                  npts_obs)
            events_tab.header['EXPSTART'] = mjd1 + delta_sec/SECPERDAY


            # DOUBLE CHECK FOR TYPO, should probably be mjd2
            delta_sec = all_delay(mjd2, parallax, objvec,
                                  earth_ephem_table,
                                  len(earth_ephem_table),
                                  obs_ephem_table,
                                  npts_obs)

            events_tab.header['EXPEND'] = mjd2 + delta_sec/SECPERDAY
            in_hdul.flush()
            if verbose:
                print("    [EVENTS,{}] extension has been updated".
                      format(e_indx))
                print("  finish time:  ", datetime.now())

        # Loop through all x1d or image extensions (there will be none if
        # the input is an events file).  Note that we only expect EXPSTART
        # and EXPEND to be present in SCI extensions; however, they
        # could be in ERR and DQ as well (e.g. in output from inttag),
        # and if they are present they must be updated.
        for e_indx in range(1, nsci_ext+1):
            cur_tab = in_hdul[e_indx]

            modified = False

            # add delaytime to EXPSTART and EXPEND, and update header
            if "EXPSTART" in cur_tab.header:
                mjd1 = cur_tab.header['EXPSTART']
                delta_sec = all_delay(mjd1, parallax, objvec,
                                      earth_ephem_table,
                                      len(earth_ephem_table),
                                      obs_ephem_table,
                                      npts_obs)
                cur_tab.header['EXPSTART'] = mjd1 + delta_sec/SECPERDAY
                modified = True

            if "EXPEND" in cur_tab.header:
                mjd2 = cur_tab.header['EXPEND']
                delta_sec = all_delay(mjd2, parallax, objvec,
                                      earth_ephem_table,
                                      len(earth_ephem_table),
                                      obs_ephem_table,
                                      npts_obs)
                cur_tab.header['EXPEND'] = mjd2 + delta_sec/SECPERDAY
                modified = True

            in_hdul.flush()
            if verbose and modified:
                print("    extension {} has been updated".format(e_indx))

        # add delaytime to TEXPSTRT and TEXPEND, and update primary header
        mjd1 = in_hdul[0].header['TEXPSTRT']
        delta_sec = all_delay(mjd1, parallax, objvec, earth_ephem_table,
                              len(earth_ephem_table), obs_ephem_table,
                              npts_obs)

        in_hdul[0].header['TEXPSTRT'] = mjd1 + delta_sec/SECPERDAY

        mjd2 = in_hdul[0].header['TEXPEND']
        delta_sec = all_delay(mjd2, parallax, objvec, earth_ephem_table,
                              len(earth_ephem_table), obs_ephem_table,
                              npts_obs)
        in_hdul[0].header['TEXPEND'] = mjd2 + delta_sec / SECPERDAY

        # add keyword to flag the fact that the times have been corrected
        in_hdul[0].header['DELAYCOR'] = ("COMPLETE",
                                         "delaytime has been applied")
        in_hdul[0].header['history'] = "Times corrected to solar system barycenter;"

        if obs_ephem_table is not None:
            in_hdul[0].header['history'] = "ORX table {}".format(obs_ephem)
        else:
            in_hdul[0].header['history'] = "no ORX tables were used."

        if verbose:
            print("... done")

        # close in file
        in_hdul.close()


def get_ephem(ephem_tables, eph_type):
    """
     read an ephemeris of state vectors
    #
    #  Description:
    #  ------------
    #  Read state vectors from one or more input tables.  The times and
    #  rectangular coordinates of the state vectors will be read from the
    #  input table from columns "TIME" (or "JD"), "X", "Y", and "Z".
    #  The times will be converted
    #  to units of Modified Julian Date, and the positions will be converted
    #  to astronomical units.
    #
    #  The number of points npts will be zero if the list of input tables
    #  is empty.  A table is required for the earth ephemeris, but none is
    #  required for the observer ephemeris.

        Note: I'm assuming the desired data is all in the
        first extension data.
    #
    #  Input table column names
    #  ------------------------
    #  JD or TIME		times corresponding to x,y,z:
    #			    in earth_ephem tables:
    #			        JD = Julian Day Number
    #			    in obs_ephem tables:
    #			        TIME = offset of time from FIRSTMJD
    #  X, Y, Z		state vector components
    #
    #  Input table parameters
    #  ----------------------
    #  "FIRSTMJD"		MJD of the first row in the table (only for obs_ephem)
    #
    #  Date		Author		Description
    #  ----		------		-----------
    #  11-Mar-1990	J.-C. Hsu	design and coding
    #  10-Jun-2003  Phil Hodge	allow multiple tables; allocate memory here;
    #				get the array of times from a table column
    #  15-Jan-2007  Phil Hodge	if keyword firstmjd is not found, read
    #				ephmstrt from primary header instead
    #------------------------------------------------------------------------------
        :param ephem_table:
        :return:
    """

    # Open the file name template string for the list of ephemeris tables.
    npts = 0
    time = None
    x, y, z = None, None, None

    if len(ephem_tables) == 0:
        return None

    # Now loop over the list of input tables and read the data.
    k = 0
    previous_mjd = 0

    for tab_file in ephem_tables:
        current_tab = Table.read(tab_file, hdu=1)
        npts += len(current_tab)

        if eph_type == EARTH_EPHEMERIS:

            # pull "JD" column
            new_time = current_tab['JD']
            new_time.name = "Time"

            # convert Julian Day Number to MJD
            new_time -= JD_TO_MJD

        else:
            # pull "TIME",
            # This isn't right I guess, can be either Time or TIME
            new_time = current_tab['Time']

            # read header parameter, the zero-point epoch
            firstmjd = g_firstmjd(tab_file)

            if firstmjd < previous_mjd:
                raise ValueError("The obs_ephem tables are out of order;\n "
                                 "use the FIRSTMJD table header keyword to "
                                 "determine the correct order.\n Times must "
                                 "be monotonically increasing")

            previous_mjd = firstmjd

            # this call pulls the units of said column into the
            # gridunit variable
            # here tooo.....
            gridunit = current_tab['Time'].unit.to_string().upper()

            # get the scale factor for converting the times to days
            if gridunit == "DAY":
                timescale = 1
            elif gridunit == "HR" or gridunit == "HOUR":
                timescale = 1 / HRPERDAY
            elif gridunit[:2] == "MIN":
                timescale = 1 / MINPERDAY
            elif gridunit[:1] == 'S':
                timescale = 1 / SECPERDAY
            else:
                raise ValueError("unrecognized TIME unit in ephemeris "
                                 "table:{}".format(gridunit))

            # scale the times, and add the zero point
            new_time = firstmjd + timescale * new_time
            new_time.name = "Time"

        #new_time.unit = cds.MJD

        if time is None:
            time = Table([new_time])
        else:
            time = vstack([time, Table([new_time])])

        x_col = current_tab['X']
        y_col = current_tab['Y']
        z_col = current_tab['Z']

        # convert the state vectors to AU
        for col in (x_col, y_col, z_col):
            colunit = col.unit.to_string().upper()
            if colunit == "AU":
                pass
            elif colunit == "KM":
                col /= KMPERAU
            else:
                raise ValueError("unrecognized {} unit in ephemeris table: "
                                 "{}".format(col.name, col.dtype))

        if x is None:
            x = Table([x_col])
            y = Table([y_col])
            z = Table([z_col])
        else:
            x = vstack([x, Table([x_col])])
            y = vstack([y, Table([y_col])])
            z = vstack([z, Table([z_col])])

        # end of loop over files

    if npts == 0:
        raise ValueError("The following input tables are empty: "
                         "{}".format(ephem_tables))

    if len(ephem_tables) > 1:
        if time['Time'][0] > time['Time'][1]:
            raise ValueError("times in {} must be monotonically "
                             "increasing.".format(ephem_tables))

    # Remove any overlapping regions.
    # loop through time column, if you hit a non-increasing value, keep
    # deleting rows until you do hit a larger number

    # Putting all columns into a new table now, so it's easier to
    # delete rows all at once.
    col_table = hstack([time, x, y, z])
    final_table = col_table[0:1]
    last_copied_time = col_table['Time'][0]

    for indx in range(1, len(col_table)):
        if last_copied_time >= col_table['Time'][indx]:
            continue
        else:
            final_table.add_row(col_table[indx])
            last_copied_time = col_table['Time'][indx]

    return final_table

def g_firstmjd(tab_file):
    '''
    g_firstmjd -- read firstmjd keyword, or use ephmstrt

    Normally the FIRSTMJD keyword will be found in the table of the ORB table.
    For data taken shortly after the ORB tables began to be written as a FITS
    table (rather than as VAX binary-format ORX files), however, FIRSTMJD was
    not yet being written into the table header.  For these files, we will
    read the EPHMSTRT keyword from the primary header and convert that string
    from SMS format to MJD.
    '''

    tab_hdulist = fits.open(tab_file)
    if 'FIRSTMJD' in tab_hdulist[1].header:
        firstmjd = tab_hdulist[1].header['FIRSTMJD']

    elif 'EPHMSTRT' in tab_hdulist[0].header:
        # Interpret ephmstrt:  year, day of year, hour, minute, second
        ephmstrt = tab_hdulist[0].header["EPHMSTRT"]
        # ephmstrt is apparently stored as a string.... so need to
        # convert to an integer
        # year.doy:hour:min:sec
        period = ephmstrt.find()
        year = int(ephmstrt[:period])
        rest = ephmstrt[period+1:].split(":")
        doy = int(rest[0])
        hour = int(rest[1])
        min = int(rest[2])
        sec = int(rest[3])
        day = doy + hour/24. + min/1400. + sec/86400.

        # Evaluate this expression with month = 1 because here 'day' is
        # the day of the year rather than day of the month.
        # reference:  Pulkinnen & VanFlandern, ApJ Suppl.
        month = 1
        k = 275 * month / 9 - 7 * (year + (month + 9) / 12) / 4 - \
            3 * ((year + (month - 9) / 7) / 100 + 1) / 4

        firstmjd = 367. * year + k + day - 678972.


    else:
        raise ValueError("Neither FIRSTMJD (ext 1) nor EPHMSTRT (ext 0) were "
                         "found in headers of file: {}".format(tab_file))

    return firstmjd


def ext_exist(ifile, in_col, nevents_tab, nsci_ext):
    """
    Check that all expected extensions are present (except GTI).
    Also check that the time column exists, if the input is a time-tag file.

    char	ifile[ARB]	# i: name of input file
    char	in_col[ARB]	# i: time column name
    int	nevents_tab	# i: number of EVENTS tables
    int	nsci_ext	# i: number of SCI (or SCI+ERR+DQ) extensions
    """

    with fits.open(ifile) as hdulist:

        # If the input is a time-tag file, look for events extensions.
        for k in range(1, nevents_tab+1):
            # can't tell exactly what i need to open, looks like
            # filename[EVENTS,k]

            if hdulist[k].header['EXTNAME'] == "EVENTS":
                if in_col not in hdulist[1].data.names:
                    print("Column {} not found in {} (skipping this file)".
                          format(in_col, ifile))
                    return False

        # If the input is a 1-D extracted spectrum or a collection of
        # image sets, open all the extensions (specified by NEXTEND).
        # for k in range(1, nsci_ext):
        # I don't think I need this...

    return True


def all_delay(epoch, parallax, objvec, earth_table, npts_earth, obs_table,
              npts_obs):
    """
    This routine computes the delaytime in seconds.

    double	epoch		# i: time (MJD) of event
    double	parallax	# i: parallax of target
    double	objvec[3]	# i: unit geocentric state vector of the target
    double	time_e[ARB], x_e[ARB], y_e[ARB], z_e[ARB]
    int	npts_earth
    double	time_o[ARB], x_o[ARB], y_o[ARB], z_o[ARB]
    int	npts_obs
    double	delta_sec	# o: correction (sec) to be added to TIME col
    #--
    char	mess[SZ_FNAME]	# for error message
    double	geomdelt, reldelt	# time corrections (sec)
    double	xyz_obs[3]	# geocentric state vector of the observer
    double	xyz_earth[3]	# baycentric state vector of the earth
    double	telvec[3]	# baycentric state vector of the observer

    """

    # calculate delaytime due to relativistic effects
    reldelt = relativ(epoch)

    # calculate the geometric delay time
    # I should really put this line into a try, add a custom intrp_state
    # Error, and catch it to give the better error message:
    # "epoch = MJD %f, outside the earth_ephem range"), -> (epoch)
    xyz_earth = intrp_state(epoch, earth_table['Time'], earth_table['X'],
                            earth_table['Y'], earth_table['Z'], npts_earth,
                            NLAN_EARTH)

    # Same goes for this guy, but with the error
    # "epoch = MJD %f, outside the obs_ephem range" -> (epoch)
    # Also I think I can take out this if zero check...
    if npts_obs > 0:
        xyz_obs = intrp_state(epoch, obs_table['Time'], obs_table['X'],
                              obs_table['Y'], obs_table['Z'], npts_obs,
                              NLAN_OBS)

    else:
        xyz_obs = [0, 0, 0]

    telvec = xyz_obs + xyz_earth
    geomdelt = geo_delay(telvec, objvec, parallax)

    delta_sec = geomdelt + reldelt

    return delta_sec


def relativ(epoch):
    """
    RELATIV -- Compute relativistic light time delay
    #
    #  Description:
    #  ------------
    #  Any clock moving with the earth is subject to an annual variation caused
    #  by two relativistic effects : (1) the time dilation of special relativity,
    #  which varies with the earth's annual variation in orbital velocity, and
    #  (2) the gravitational redshift of general relativity, according to which
    #  a clock runs more slowly in a stronger gravitational field.
    #
    #  RELATIV performs the transformation from proper time on earth to coordinate
    #  time in space-time frame of reference in which the solar system barycenter
    #  is at rest.  The algorithm is from Moyer, Cel. Mech., 23, 33, 1981, but
    #  with the changing value of the Earth orbit eccentricity taken into account
    #  together with a more accurate value for the coefficient of the main annual
    #  term, and with two extra terms included to allow for the eccentricities of
    #  the orbits of Jupiter and Saturn.  It is known that Moyer's equation is
    #  accurate to about a 7 microsecond level.
    #
    #  In practice, to better than sufficient accuracy it is enough to use the
    #  proper time on earth (ephemeris time, terrestrial time (TT) or atomic
    #  time (TAI); the small offset of 32.184 seconds doesn't matter.)
    #  The output time difference
    #
    #	DELTA_T = T - TAI ,
    #
    #  where T is the coordinate time in the solar system barycenter frame of
    #  reference, and TAI is the International Atomic Time from the clock on
    #  earth.  DELTA_T has a main (annual) term of amplitude approximately
    #  1.7 milliseconds.
    #
    #  References:
    #  Escobal, P.  1967, Methods of Astrodynamics, Wiley, p. 8.
    #  Moyer, T.    1981, Celestial Mechanics, vol. 23, pp. 33 & 57.
    #
    #  Date		Author			Description
    #  ----		------			-----------
    #  11-Dec-1984	C. D. Biemesderfer	Modified Starlink RCC2 function
    #  11-Mar-1990	J.-C. Hsu		rewrite in SPP
    #   6-Sep-2005	Phil Hodge		Remove "/ RADIAN" factors;
    #					don't include <math.h>

    double	epoch		# input: MJD of the observation
    double	delta_t		# output: Relativistic time delay (sec)

    # Assorted orbital parameters
    double	ecc_e		# Eccentricity of Earth orbit
    double	e_anom_e	# Eccentric anomaly of Earth
    double	m_anom_e	# Mean anomaly of Earth
    double	m_elon_m	# Mean elongation of Moon from Sun

    double	m_anom_j	# Mean anomaly of Jupiter
    double	l_m_lj		# Mean longitude of Jupiter from E-M barycenter
    double	omega_j		# Related to argument of perihelion for Jupiter

    double	m_anom_s	# Mean anomaly of Saturn
    double	l_m_ls		# Mean longitude of Saturn  from E-M barycenter
    double	omega_s		# Related to argument of perihelion for Saturn

    double	time		# Obs time in secs past 1950 Jan 1.0
    """

    # Convert input time to seconds past 1950 Jan 1.0
    time = (epoch - 33282.) * 86400.

    # Calculate orbital parameters at time of observation
    ecc_e = 1.673014e-2 - 1.325e-14 * time        # Escobal (1.5)

    m_anom_e = 6.248291  + 1.99096871e-7 * time    # Moyer II (44)
    m_elon_m = 2.518411  + 2.462600818e-6 * time	# Moyer II (45)
    l_m_lj = 5.652593  + 1.82313637e-7 * time	    # Moyer II (46)
    l_m_ls = 2.125474  + 1.92339923e-7 * time	    # Moyer II (47)
    m_anom_j = 5.286877  + 1.6785063e-8 * time	    # Moyer II (48)
    m_anom_s = 1.165341  + 0.6758558e-8 * time	    # Moyer II (49)

    e_anom_e = m_anom_e + ecc_e * sin(m_anom_e)     # Moyer II (40)

    omega_j = l_m_lj - m_anom_j  # Added by Starlink
    omega_s = l_m_ls - m_anom_s  # Added by Starlink

    # Correction from proper time to coordinate time, using Moyer II (38)

    # Daily motion of Earth around Sun correction has amplitude of millisecs

    delta_t = (9.9153e-2 * ecc_e) * sin(e_anom_e) + \
               1.548e-6 * sin (m_elon_m) + \
               5.21e-6 * sin(m_anom_j) + \
     		   2.45e-6 * sin(m_anom_s) + \
     		   20.73e-6 * sin(l_m_lj) + \
     		   4.58e-6 * sin(l_m_ls) + \
     		   1.00e-6 * sin(omega_j) + 0.26e-6 * sin(omega_s) # Added by Starlink

    return delta_t


def intrp_state(epoch, time, x, y, z, npts, nlag):
    """
    #  INTRP_STATE -- state vector interpolation for a specified epoch
    #
    #  Description:
    #  ------------
    #  Using Lagrange's formula of polynomial interpolation, obtain the state
    #  vector of any epoch from the input state vector table.
    #  The input array of times must be monotonically increasing, but they
    #  don't have to be uniformly spaced.
    #
    #  Date 	Author 		Description
    #  ----		------		-----------
    #  28-Feb-1988  J.-C. Hsu       original, vxyzin.for
    #  17-Apr-1990	J.-C. Hsu	rewrite in SPP
    #  10-Jun-2003  Phil Hodge	include time in calling sequence, instead of
    #				computing time from first value and increment

    double	epoch		# i: desired epoch
    double	time[npts]	# i: times (MJD) corresponding to x, y, z
    double	x[npts], y[npts], z[npts]	# i: state vectors
    double	xyz[3]		# o: state vector of the desired epoch
    int	npts		# i: number of entries of the input state vectors
    int	nlag		# i: number of points used in Lagrange's
                #        polynomial interpolation formula

    int	icenter, istart, istop
    int	i1, i2		# indexes for searching for epoch in time array
    double	dxyz[3]		# error of xyz (ignored)
    pointer	sp, c, d	# work space for vpolin
    """

    # find the array index of the desired epoch (binary search)
    i1 = 1
    i2 = npts

    while (i2 - i1) > 1:
        icenter = int((i1 + i2) / 2)
        if epoch > time[icenter]:
            i1 = icenter
        else:
            i2 = icenter

    icenter = i1

    istart = icenter - int(nlag / 2)
    istop = istart + nlag - 1

    # make sure the starting and stopping points are within limits
    if istart < 1 or istop > npts:
        raise ValueError("epoch = {:.4f}; ephemeris range = {:.4f} to {:.4f} \n"
                         "Input epoch is out of ephemeris range".
                         format(epoch, time[0], time[npts-1]))

    # perform the interpolation
    # time (x input)  x (y input)  nlag (size of arrays)
    # epoch (x value to be calc) xyz[1] (output y) dxyz[1] (error)
    xyz = np.zeros((3))
    time_in = time[istart-1: istart+nlag]
    xyz[0] = vpolin(time_in, x[istart-1: istart+nlag], epoch)
    xyz[1] = vpolin(time_in, y[istart-1: istart+nlag], epoch)
    xyz[2] = vpolin(time_in, z[istart-1: istart+nlag], epoch)

    return xyz


def vpolin(xa, ya, in_x):
    """Lagrange interpolation on input data arrays"""

    # try and use scipy lagrange interpolation here,
    # looks like it is a implementation of neville's algorithm
    # It's taking nlag=10 numbers, assuming it's pulling those from
    # the start of the array.... although this is a bit dangerous
    poly = interp1d(xa, ya, kind="cubic")
    #poly = lagrange(xa, ya)

    return poly(in_x)


def geo_delay(stvec, objvec, parallax):
    """
    #  GEO_DELAY -- Get time delay due to displacement from solar-system barycenter
    #
    #  Description:
    #  ------------
    #  Find the difference between actual time of observation and the time that
    #  would have been observed at the solar-system barycenter.  The object which
    #  was observed is assumed to be very distant (outside the solar system).
    #
    #  This module applies the geometric correction due to displacement, not
    #  the correction due to gravitational redshift.  Geometric effects include
    #  the corrections for Earth orbital motion (amplitude of 499 sec) and for
    #  ST orbital motion (amplitude of 23 millisec)
    #
    #  A comment on accuracy required:  one arcsec difference in the position of
    #  the object can make a difference in timing of 2.5 millisec.
    #
    #  Date		Author			Description
    #  ----		------			-----------
    #  15-Nov-1984	C. D. Biemesderfer	Original module
    #  11-Apr-1990	J.-C. Hsu		rewrite in SPP

    double	stvec[3]	# input: telescope state vector (AU)
    double	objvec[3]	# input: unit state vector of the target (unitless)
    double	parallax	# input: Trigonometric parallax (")
    double	geomdelt	# output: Geometric time correction (sec)

    double	ob_dot_st	# inner product of OBJVEC and STVEC
    double	st_dot_st	# Squared modulus of STVEC
    double	pxau		# Parallax expressed in AU
    double	rsqrterm	# 2nd order term in (r_cos_theta - amplitude of 499 s)

    double	adotd()
    """

    # Compute inner product of object unit state vector and telescope state
    # vector
    ob_dot_st = np.dot(objvec, stvec)    # in AU
    st_dot_st = np.dot(stvec, stvec)     # in AU**2

    # Compute parallax over 2 coefficient in AU and term involving
    # the square of the telescope's barycentric distance
    pxau = parallax / AUPERPC
    rsqrterm = st_dot_st - ob_dot_st * ob_dot_st

    # Low precision correction includes first two terms in series expansion.
    # Amplitude of second term is on the order of 1 millisec.
    geomdelt = ob_dot_st - (pxau / 2.0) * rsqrterm

    # This was commented out in the code during time of port
    """
    #if (precision .eq. high)

	# Include third term (third order in r_cosine_theta - amplitude of order
	# nanosec or smaller for extra-solar-system objects although much larger
	# for solar system observation)
	    #geomdelt = geomdelt - (pxau * pxau / 2.d0) * ob_dot_st * rsqrterm
    """

    # Correct result in AU to sec
    geomdelt = geomdelt * CLIGHT

    return geomdelt
