#! /usr/bin/env python

from astropy.io import fits
from astropy.table import Table

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
KMPERPC	= 3.085678   #D13	# kilometers per parsec
AUPERPC	= 206264.8062470964  #D0	# how many AU in one parsec
KMPERAU	= KMPERPC/AUPERPC
CLIGHT = 499.00479   #D0	# light travel time (in sec) of 1 AU

# might not need these
EARTH_EPHEMERIS = 1
OBS_EPHEMERIS = 2

JD_TO_MJD = 2400000.5   #d0	# subtract from JD to get MJD



def odelaytime(table_names, earth_ephem, obs_ephem, distance, dist_unit,
               in_col, verbose=False):

    parallax = distance

    # convert distance to parsec first
    if dist_unit == "au":
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
        raise ValueError("non-postive distance: {}".format(parallax))

    # I don't think I need this part
    # call smark(sp)
    # call salloc(history, SZ_FNAME, TY_CHAR)
    # call salloc(extname, SZ_FNAME, TY_CHAR)





# odelay_do
# call odelay_do (fin, nfin, parallax, earth_ephem, obs_ephem,
# in_col, verbose)

def get_ephem(ephem_tables, eph_type):
    '''
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
    '''

    # Open the file name template string for the list of ephemeris tables.
    npts = 0
    time = 0
    x, y, z = 0, 0, 0

    if len(ephem_tables) == 0:
        return npts, time, x, y, z

    # Count the total number of rows in all the input tables.
    # Not sure if I need this
    #nrows = tbpsta (tp, TBL_NROWS)
    #npts = npts + nrows

    if npts == 0:
        raise ValueError("The following input tables are empty: "
                         "{}".format(ephem_tables))

    # Now loop over the list of input tables and read the data.
    k = 0
    previous_mjd = 0

    for tab_file in ephem_tables:
        # open table file
        #tp = tbtopn (table, READ_ONLY, 0)
        current_tab = Table.read(tab_file, hdu=1)

        if eph_type == EARTH_EPHEMERIS:

            #pull "JD" column?
            # call tbcfnd1 (tp, "JD", cp[1])
            new_time = current_tab['JD']

            # convert Julian Day Number to MJD
            new_time -= JD_TO_MJD

        else:
            #pull "TIME",
            # call tbcfnd1 (tp, "TIME", cp[1])
            new_time = current_tab['TIME']

            # scale the times, and add the zero point
            new_time = firstmjd + timescale * new_time


        if time == 0:
            time = new_time
        else:
            time.insert(-1, new_time)

        # pull "X", "Y", and "Z" columns
	    #call tbcfnd1 (tp, "X", cp[2])
	    #call tbcfnd1 (tp, "Y", cp[3])
	    #call tbcfnd1 (tp, "Z", cp[4])
        if x == 0:
            x = current_tab['X']
            y = current_tab['X']
            z = current_tab['Z']
        else:
            x.insert(-1, current_tab['X'])
            y.insert(-1, current_tab['Y'])
            z.insert(-1, current_tab['Z'])

	    if eph_type == OBS_EPHEMERIS:

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
            #call tbcigt (cp[1], TBL_COL_UNITS, gridunit, SZ_COLUNITS)
            gridunit = current_tab['TIME'].unit.to_string().upper()

            # get the scale factor for converting the times to days
            if gridunit == "DAY":
                timescale = 1
            elif gridunit == "HR" or gridunit == "HOUR":
                timescale = 1 / HRPERDAY
            elif gridunit[:2] == "MIN":
                timescale = 1 / MINPERDAY
            elif gridunit[1] == 'S':
                timescale = 1 / SECPERDAY
            else:
                raise ValueError("unrecognized TIME unit in ephemeris "
                                 "table:{}".format(gridunit))


	    # read the data (already did this)
        # tbcgtd - Read values for one column from a range of rows.
	    #  time (or mjd) = cp[1] = Memd[time+k]
	    #  x = cp[2] = Memd[x+k]
	    #  y = cp[3] = Memd[y+k]
	    #  z = cp[4] = Memd[z+k]

        # convert the state vectors to AU
        for col in (x,y,z):
            colunit = col.unit.to_string().upper()
            if colunit == "AU":
                pass
            elif colunit == "KM":
                #call adivkd (Memd[x+k], KMPERAU, Memd[x+k], nrows)
                # ADIVK -- Divide a vector by a constant
                col = col/KMPERAU
            else:
                raise ValueError("unrecognized {} unit in ephemeris table: "
                                 "{}".format(col.name, col.dtype))

        # end of loop over files

    if len(ephem_tables) > 1:
        if [time][0] > time[1]:
            raise ValueError("times in {} must be monotonically "
                             "increasing.".format(ephem_tables))

    # Remove any overlapping regions.
    # loop through time column, if you hit a non-increasing value, keep
    # deleting rows until you do hit a larger number

    # Putting all columns into a new table now, so it's easier to
    # delete rows all at once.
    final_table = Table()




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
