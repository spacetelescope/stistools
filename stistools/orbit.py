import math
from astropy.io import fits

TWOPI = (math.pi * 2.0)
SEC_PER_DAY = 86400.0


class HSTOrbit(object):
    """Orbital parameters.

    The public methods are getOrbitper and getPos.
    """

    def __init__(self, spt):
        """Orbital parameters.

        Parameters
        ----------
        spt: str
            Name of the support file (rootname_spt.fits).
        """

        self.orb = {}
        self._readOrbitalParameters(spt)

    def _readOrbitalParameters(self, spt):
        """Get the orbital parameters from the spt primary header.

        Parameters
        ----------
        spt: str
            Name of the support file (rootname_spt.fits)
        """

        fd = fits.open(spt, mode="readonly")
        phdr = fd[0].header

        self.orb = {}
        self.orb["argperig"] = phdr["argperig"]
        self.orb["cirveloc"] = phdr["cirveloc"]
        self.orb["cosincli"] = phdr["cosincli"]
        self.orb["ecbdx3"] = phdr["ecbdx3"]
        self.orb["eccentry"] = phdr["eccentry"]
        self.orb["eccentx2"] = phdr["eccentx2"]
        self.orb["ecbdx4d3"] = phdr["ecbdx4d3"]
        self.orb["epchtime"] = phdr["epchtime"]
        self.orb["esqdx5d2"] = phdr["esqdx5d2"]
        self.orb["fdmeanan"] = phdr["fdmeanan"]
        self.orb["hsthorb"] = phdr["hsthorb"]
        self.orb["meananom"] = phdr["meananom"]
        self.orb["rascascn"] = phdr["rascascn"]
        self.orb["rcargper"] = phdr["rcargper"]
        self.orb["rcascnrv"] = phdr["rcascnrv"]
        self.orb["sdmeanan"] = phdr["sdmeanan"]
        self.orb["semilrec"] = phdr["semilrec"]
        self.orb["sineincl"] = phdr["sineincl"]

        fd.close()

    def getOrbitper(self):
        """Return the orbital period.

        Returns
        -------
        orbit period: float
            Orbital period in seconds.
        """

        return 2. * self.orb["hsthorb"]

    def getPos(self, mjd):
        """Get position and velocity at a given time.

        # S. Hulbert, Oct 91    Original
        # PEH, 2008 Oct 3       Converted from SPP to Python

        Parameters
        ----------
        mjd: float
            Time (Modified Julian Date).

        Returns
        -------
        (x_hst, v_hst): tuple of two vectors (3-element lists)
            Position and velocity at the specified time.
        """

        # These will be returned, after assigning the actual values.
        x_hst = [0., 0., 0.]
        v_hst = [0., 0., 0.]

        argperig = self.orb["argperig"]
        cirveloc = self.orb["cirveloc"]
        cosincli = self.orb["cosincli"]
        ecbdx3 = self.orb["ecbdx3"]
        eccentry = self.orb["eccentry"]
        eccentx2 = self.orb["eccentx2"]
        ecbdx4d3 = self.orb["ecbdx4d3"]
        epchtime = self.orb["epchtime"]
        esqdx5d2 = self.orb["esqdx5d2"]
        fdmeanan = self.orb["fdmeanan"]
        meananom = self.orb["meananom"]
        rascascn = self.orb["rascascn"]
        rcargper = self.orb["rcargper"]
        rcascnrv = self.orb["rcascnrv"]
        sdmeanan = self.orb["sdmeanan"]
        semilrec = self.orb["semilrec"]
        sineincl = self.orb["sineincl"]

        # convert time from MJD to seconds since 1985 Jan 1
        sec85 = (mjd - 46066.0) * SEC_PER_DAY

        # calculate time difference between observation and epoch time
        deltim = sec85 - epchtime

        # mean anomaly
        temp2 = fdmeanan * deltim
        temp3 = 0.5 * sdmeanan * deltim*deltim
        m = meananom + TWOPI * (temp2 + temp3)

        sin_m = math.sin(m)
        cos_m = math.cos(m)

        # true anomaly (equation of the center)
        v = m + sin_m * (eccentx2 + ecbdx3 * cos_m * cos_m -
                         ecbdx4d3 * sin_m * sin_m + esqdx5d2 * cos_m)
        sin_v = math.sin(v)
        cos_v = math.cos(v)

        # distance
        r = semilrec / (1.0 + eccentry * cos_v)

        # argument of perigee
        wsmall = TWOPI * (argperig + rcargper * deltim)

        # longitude of the ascending node
        wbig = TWOPI * (rascascn + rcascnrv * deltim)
        sin_wbig = math.sin(wbig)
        cos_wbig = math.cos(wbig)

        # calculate the rectangular coordinates
        #  (see Smart, Spherical Astronomy, section 75, page 122-124)

        f = wsmall + v
        sin_f = math.sin(f)
        cos_f = math.cos(f)

        x_hst[0] = r * (cos_wbig * cos_f - cosincli * sin_wbig * sin_f)
        x_hst[1] = r * (sin_wbig * cos_f + cosincli * cos_wbig * sin_f)
        x_hst[2] = r * sineincl * sin_f

        a0 = cirveloc * eccentry * sin_v / r
        a1 = cirveloc * (1.0 + eccentry * cos_v) + \
            TWOPI * rcargper * r
        v_hst[0] = a0 * x_hst[0] - \
            a1 * (cos_wbig * sin_f + cosincli * sin_wbig * cos_f) - \
            TWOPI * rcascnrv * x_hst[1]
        v_hst[1] = a0 * x_hst[1] - \
            a1 * (sin_wbig * sin_f - cosincli * cos_wbig * cos_f) + \
            TWOPI * rcascnrv * x_hst[0]
        v_hst[2] = a0 * x_hst[2] + a1 * sineincl * cos_f

        # Convert from meters to kilometers.
        for i in range(3):
            x_hst[i] /= 1000.0
            v_hst[i] /= 1000.0

        return x_hst, v_hst
