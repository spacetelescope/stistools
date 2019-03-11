#! /usr/bin/env python

import sys
import math
import numpy as np

from astropy.io import fits

from . import observation
from . import orbit

__doc__ = """
This class computes Doppler shift information for each imset of a dataset.
Keywords will be read from the science file and from the support file. The
Doppler shift information is printed to the standard output by default, but
this can be turned off by setting quiet to True. Three task parameters will be
updated by doppinfo; these allow the user to compute the Doppler shift in
high-res pixels or km/s at any time close to the time of the exposure.

The printed information will be in one of two formats, depending on the value
of increment. If increment is zero, the average and extreme values of the
Doppler shift during the exposure will be printed. If increment is
greater than zero, the actual Doppler shift and radial velocity will be
printed at the beginning of the exposure and at every increment seconds
thereafter until the end of the exposure. Both the printed value of the
Doppler shift and the doppmag parameter will be in high-res pixels.

Most of the Doppler shift information is computed directly from the orbital
elements of HST, as given in the support file primary header. Some information,
however, is computed based on the approximation of a circular orbit. This
approximation is used for the average Doppler shift during the exposure
(printed if increment is zero) and for the task parameters doppzero, doppmag
and radvel that are updated by doppinfo. These parameters are applied as
terms in a sine function, which inherently involves a circular-orbit
approximation.

The parameters for the circular orbit are determined as follows. The HST
orbital elements are gotten from the primary header of the support file. The
target position is taken from the keywords RA_TARG and DEC_TARG in the science
file primary header. The velocity of HST is computed from the orbital elements
at 64 equally spaced times throughout an orbit, centered on the midpoint of the
exposure, (EXPSTART + EXPEND) / 2, and the component of this velocity in the
direction away from the target (i.e. the radial velocity) is taken. A sine
function is fit to these radial velocities; the amplitude is radvel, and the
amplitude and phase are used to compute doppmag and doppzero.


Examples
--------

:class:`Doppinfo` with dt of 100:

>>> import stistools
>>> stistools.doppinfo.Doppinfo("ocb6o2020_raw.fits", dt=100, spt="ocb6o2020_spt.fits")
# orbitper  doppzero      doppmag      doppmag_v    file
  5728.67  56752.114170  11.68643135   7.40391177   ocb6o2020_raw.fits[sci,1]
# time (MJD)   shift   radvel
56752.165175  -11.59   -7.345
56752.166333  -11.37   -7.203
56752.167490  -11.01   -6.975
56752.168647  -10.52   -6.663
56752.169805   -9.90   -6.272
56752.170962   -9.16   -5.805
56752.172120   -8.32   -5.269
56752.173277   -7.37   -4.669
56752.174434   -6.34   -4.014
56752.175592   -5.23   -3.311
56752.176749   -4.05   -2.568
56752.177907   -2.83   -1.794
56752.179064   -1.58   -0.998
56752.180222   -0.30   -0.190
# orbitper  doppzero      doppmag      doppmag_v    file
  5728.67  56752.180505  11.68734454   7.40449032   ocb6o2020_raw.fits[sci,2]
# time (MJD)   shift   radvel
56752.181784    1.42    0.902
56752.182941    2.68    1.700
56752.184099    3.91    2.477
56752.185256    5.09    3.225
56752.186413    6.21    3.935
56752.187571    7.26    4.598
56752.188728    8.22    5.205
56752.189886    9.08    5.750
56752.191043    9.83    6.227
56752.192200   10.46    6.628
56752.193358   10.97    6.950
56752.194515   11.35    7.189
56752.195673   11.59    7.342
56752.196830   11.69    7.406


"""

__version__ = "3.0"
__vdate__ = "2018-09-14"

# multiply by DEG_RAD to convert from degrees to radians
DEG_RAD = (math.pi / 180.0)
TWOPI = (math.pi * 2.0)
SPEED_OF_LIGHT = 299792.458             # km / s
SEC_PER_DAY = 86400.0

# number of points for computing circular-orbit parameters
NPTS = 64
NPTS_D = 64.


class Doppinfo(object):
    """Compute Doppler parameters and information from HST orbital elements.
    This class previously supported both COS and STIS data, but now only
    supports STIS data.  The class will print doppler shift information for all
    imsets contained in the input image.

    Results will be printed to standard out.  To have the DOPPZERO, DOPPMAG,
    and DOPPPMAGV keywords inserted/updated in the header you can set update to
    True."""

    def __init__(self, input, spt=None, dt=0., update=False, quiet=False):
        """Compute Doppler parameters.

        Parameters
        ----------
        input: str
            The name of an input file; Only STIS data is supported.

        spt: str or None
            The name of the support (_spt.fits) file to use for getting
            the orbital parameters; if not specified, the spt file name
            will be constructed from the input file name by truncating
            after "raw", "corrtag", "flt", "counts" or "x1d" and appending
            "spt.fits".

        dt: float
            Time interval (seconds) for printing Doppler shift throughout
            the orbit.

        update: boolean
            By default, Doppinfo just prints the Doppler parameter values.
            If you specify update=True, the input files will be modified
            (in-place) by updating the keywords ORBITPER, DOPPZERO,
            DOPPMAG, and DOPPMAGV in the first extension.

        quiet: boolean
            If True, values of Doppler parameters and other info will not
            be printed.
        """

        self.input = input
        self.update = update
        self.quiet = quiet
        self.obs = {}

        # unit vector (3-element list) pointing toward the target
        self.target = None

        instrument, nextend = self._getInitInfo()
        self.sci_num = nextend // 3
        dt /= SEC_PER_DAY  # convert dt to days

        for sci_ext in np.arange(self.sci_num)+1:

            self.obs = observation.initObservation(input, instrument, sci_ext)
            self.obs.getInfo()

            if spt is None:
                spt = self._findSptName()

            self.orbit = orbit.HSTOrbit(spt)

            self._getDoppParam()

            if not quiet:
                print("# orbitper  doppzero      doppmag      doppmag_v    "
                      "file")
                print("  {:.7g}  {:.6f}  {:.8f}   {:.8f}   {}[sci,{}]".
                      format(self.orbitper, self.doppzero, self.doppmag,
                             self.doppmag_v, self.input.split("/")[-1],
                             sci_ext))
                self.printDopplerShift(dt)

            # Maybe shouldn't be in loop
            if update:
                self._updateKeywords(input, sci_ext)

    def _getInitInfo(self):
        """
        Get nextend and instrument.
        """

        fd = fits.open(self.input, mode="readonly")
        instrument = fd[0].header.get("instrume", "missing")
        nextend = fd[0].header['nextend']
        fd.close()

        return instrument, nextend

    def _findSptName(self):
        """Get the name of the support file.

        Returns
        -------
        spt: str
            Name of the support (_spt.fits) file.
        """

        i = self.input.rfind("raw")
        if i < 0:
            i = self.input.rfind("corrtag")
            if i < 0:
                i = self.input.rfind("flt")
                if i < 0:
                    i = self.input.rfind("counts")
                    if i < 0:
                        i = self.input.rfind("x1d")
        if i >= 0:
            spt = self.input[0:i] + "spt.fits"
        else:
            raise RuntimeError("Don't understand input file name '{}'".
                               format(self.input))

        return spt

    def _getDoppParam(self):
        """Compute Doppler parameters.

        The following attributes will be assigned:
        orbitper:  orbital period (seconds)
        doppzero:  time (MJD) when the Doppler shift is zero and increasing
        doppmag_v:  magnitude of the Doppler shift in km/s
        doppmag:  magnitude of the Doppler shift in pixels
        """

        # Convert target ra,dec to rectangular coordinates (unit radius).
        self.target = self._sph_rec(self.obs.ra_targ * DEG_RAD,
                                   self.obs.dec_targ * DEG_RAD)

        self.orbitper = self.orbit.getOrbitper()         # seconds
        orbit_period = self.orbitper / SEC_PER_DAY       # days

        expmiddle = (self.obs.expstart + self.obs.expend) / 2.0

        # Compute Fourier coefficients based on NPTS points in an orbit,
        # centered on the middle of the exposure.

        sum_sin = 0.0
        sum_cos = 0.0
        t_origin = expmiddle - orbit_period / 2.0
        for i in range(NPTS):
            delt = i * (orbit_period / NPTS_D)
            time = t_origin + delt
            radvel = self._get_rv(time)
            sum_sin += radvel * math.sin(TWOPI * delt / orbit_period)
            sum_cos += radvel * math.cos(TWOPI * delt / orbit_period)

        # Normalize by dividing by the sum of (sin**2) at NPTS equally spaced
        # times in one orbit.

        acoeff = sum_sin / (NPTS_D/2.0)
        bcoeff = sum_cos / (NPTS_D/2.0)

        """Find doppzero and doppmag, assuming a circular orbit.

         Assume that the radial velocity has this form:
           radvel = acoeff * sin(arg) + bcoeff * cos(arg)
         where:
           arg = (time - t_origin) * 2*pi/P
           time is in MJD
           t_origin = MJD at the middle of the exposure - P/2
           P is the orbital periond of HST
         Write:
           acoeff = doppmag_v * cos(theta)
           bcoeff = doppmag_v * sin(theta)
         theta will give us doppzero, as explained below.
         Then the expression for radial velocity vs time is:
           radvel = doppmag_v * sin(arg + theta)
         From the definitions of doppmag_v and doppzero:
           radvel = doppmag_v * sin((time - doppzero) * 2*pi/P)
         so:
           arg + theta = (time - doppzero) * 2*pi/P
           (time - t_origin) * 2*pi/P + theta = (time - doppzero) * 2*pi/P
           theta = [(time - doppzero) - (time - t_origin)] * 2*pi/P
                 = (t_origin - doppzero) * 2*pi/P
         and:
           doppzero = -theta * P/(2*pi) + t_origin
         We have already computed acoeff and bcoeff, so:
           theta = atan2(bcoeff, acoeff)
         """

        self.doppzero = -math.atan2(bcoeff, acoeff) \
            * orbit_period / TWOPI + t_origin
        self.doppmag_v = math.sqrt(acoeff*acoeff + bcoeff*bcoeff)
        self.doppmag = self._rvToPixels(self.doppmag_v)

    def _sph_rec(self, ra_targ, dec_targ):
        """Convert from RA & Dec to rectangular coordinates.

        Parameters
        ----------
        ra_targ: float
            Right ascension (radians) of the target.

        dec_targ: float
            Declination (radians) of the target.

        Returns
        -------
        target: list of three floats
            Unit vector (rectangular coords.) pointing toward the target.
        """

        target = [0., 0., 0.]
        target[0] = math.cos(dec_targ) * math.cos(ra_targ)
        target[1] = math.cos(dec_targ) * math.sin(ra_targ)
        target[2] = math.sin(dec_targ)

        return target

    def _get_rv(self, time):
        """Compute the radial velocity.

        Parameters
        ----------
        time: float
            A particular time (MJD).

        Returns
        -------
        radial velocity: float
            The component of the velocity away from the target.
        """

        (x, v) = self.orbit.getPos(time)

        # This is the component of velocity toward the target.
        dot_product = self.target[0] * v[0] + self.target[1] * v[1] + \
            self.target[2] * v[2]

        # Change the sign to get the component away from the target.
        return -dot_product

    def _rvToPixels(self, radvel):
        """Convert radial velocity to Doppler shift in pixels.

        Parameters
        ----------
        radvel: float
            Radial velocity in km/s.

        Returns
        -------
        doppmag: str
            Maximum value of Doppler shift (pixels) over entire orbit.
        """

        doppmag = radvel * self.obs.cenwave / \
            (SPEED_OF_LIGHT * self.obs.dispersion)
        return doppmag

    def _pixelsToRv(self, doppmag):
        """Convert Doppler shift in pixels to radial velocity.

        Parameters
        ----------
        doppmag: float
            Maximum value of Doppler shift (pixels) over entire orbit.

        Returns
        -------
        radvel: float
            Radial velocity in km/s.
        """

        radvel = doppmag * SPEED_OF_LIGHT * self.obs.dispersion / \
            self.obs.cenwave
        return radvel

    def printDopplerShift(self, dt):
        """Compute and print the Doppler shift at intervals of dt.

        Parameters
        ----------
        dt: float
            Time interval (seconds) for printing Doppler shift throughout
            the orbit, or if dt is zero print the min and max Doppler shift
            during the orbit.
        """

        expstart = self.obs.expstart
        expend = self.obs.expend
        expmiddle = (expstart + expend) / 2.
        doppzero = self.doppzero
        orbit_period = self.orbitper / SEC_PER_DAY       # days

        if dt > 0.:
            print("# time (MJD)   shift   radvel")

            # Add 1.e-4 to expend to include end of interval, in case
            # increment divides exposure time evenly.
            done = False
            time = expstart
            while not done:
                if time <= expend+1.e-4:
                    radvel = self._get_rv(time)
                    doppmag = self._rvToPixels(radvel)
                    print("{:12.6f} {:7.2f} {:8.3f}".
                          format(time, doppmag, radvel))
                    time += dt
                else:
                    done = True
        else:
            # Use the radial velocity at the middle of the exposure
            # as an initial value for finding min & max radial velocity.
            mid_radvel = self._get_rv(expmiddle)
            min_radvel = mid_radvel
            max_radvel = mid_radvel
            t_min = expmiddle
            t_max = expmiddle
            delta = 1.e-4
            done = False
            time = expstart + delta
            while not done:
                if time <= expend:
                    radvel = self._get_rv(time)
                    if radvel < min_radvel:
                        min_radvel = radvel
                        t_min = time
                    if radvel > max_radvel:
                        max_radvel = radvel
                        t_max = time
                    time += delta
                else:
                    done = True
            # Explicitly check the radial velocities at the endpoints.
            min_at_end = False                        # initial values
            max_at_end = False
            radvel = self._get_rv(expstart)
            if radvel < min_radvel:
                min_radvel = radvel
                t_min = expstart
                min_at_end = True
            if radvel > max_radvel:
                max_radvel = radvel
                t_max = expstart
                max_at_end = True
            radvel = self._get_rv(expend)
            if radvel < min_radvel:
                min_radvel = radvel
                t_min = expend
                min_at_end = True
            if radvel > max_radvel:
                max_radvel = radvel
                t_max = expend
                max_at_end = True
            # Improve the values of min and max radial velocity.
            rv = [0., 0., 0.]
            if not min_at_end:
                rv[0] = self._get_rv(t_min-delta)
                rv[1] = self._get_rv(t_min)
                rv[2] = self._get_rv(t_min+delta)
                time = self._peakQuadratic(rv, t_min, delta)
                time = max(time, expstart)
                time = min(time, expend)
                min_radvel = self._get_rv(time)

            if not max_at_end:
                rv[0] = self._get_rv(t_max-delta)
                rv[1] = self._get_rv(t_max)
                rv[2] = self._get_rv(t_max+delta)
                time = self._peakQuadratic(rv, t_max, delta)
                time = max(time, expstart)
                time = min(time, expend)
                max_radvel = self._get_rv(time)

            # Compute the average radial velocity.  Note that this
            # assumes a circular orbit.
            if expend == expstart:
                avg_radvel = self._get_rv(expmiddle)
            else:
                avg_dopp = self.doppmag * \
                    (math.cos(TWOPI * (expstart - doppzero) / orbit_period) -
                    math.cos(TWOPI * (expend - doppzero) / orbit_period)) * \
                    orbit_period / TWOPI / (expend - expstart)
                avg_radvel = self._pixelsToRv(avg_dopp)
            mid_dopp = self._rvToPixels(mid_radvel)
            avg_dopp = self._rvToPixels(avg_radvel)
            min_dopp = self._rvToPixels(min_radvel)
            max_dopp = self._rvToPixels(max_radvel)
            print("# midpoint   midpoint Doppler  average Doppler  "
                  "minimum Doppler  maximum Doppler")
            print("#   MJD        pixels   km/s    pixels   km/s    "
                  "pixels   km/s    pixels   km/s")
            print("{:12.6f} {:8.2f} {:6.3f}  {:8.2f} {:6.3f}  "
                  "{:8.2f} {:6.3f}  {:8.2f} {:6.3f}  {}".
                  format(expmiddle, mid_dopp, mid_radvel, avg_dopp, avg_radvel,
                         min_dopp, min_radvel, max_dopp, max_radvel,
                         self.input.split("/")[-1]))

        print("")

    def _peakQuadratic(self, y, x_middle, spacing):
        """Get the location of the maximum (or minimum) of a quadratic.

        Parameters
        ----------
        y: array_like
            Values of a function at three uniformly spaced points.

        x_middle: float
            Independent variable at the middle point.

        spacing: float
            Increment in the independent variable between elements of `y`.

        Returns
        -------
        x: float
            Independent variable of the maximum (or minimum) of the
            quadratic that passes through the three uniformly spaced
            points.
        """

        denominator = y[0] - 2.0 * y[1] + y[2]
        if denominator == 0.0:
            return x_middle

        dx = (y[0] - y[2]) / (2.0 * denominator)

        return dx * spacing + x_middle

    def _updateKeywords(self, input, sci_ext):
        """Update keywords in the first extension header.

        Parameters
        ----------
        input: str
            The name of an input file (modified in-place).

        sci_ext: int
            The number of the science extension.
        """

        fd = fits.open(input, mode="update")
        hdr = fd['sci', sci_ext].header

        old_orbitper = hdr.get("orbitper", -999)
        old_doppzero = hdr.get("doppzero", -999)
        old_doppmag = hdr.get("doppmag", -999)
        old_doppmag_v = hdr.get("doppmagv", -999)

        hdr["orbitper"] = self.orbitper
        hdr["doppzero"] = self.doppzero
        hdr["doppmag"] = self.doppmag
        hdr["doppmagv"] = self.doppmag_v

        fd.close()

        if not self.quiet:
            print("{}[sci,{}] has been updated as follows:".
                  format(input.split("/")[-1], sci_ext))
            if old_orbitper == -999:
                print("orbitper:  {:.4f} (added)".format(self.orbitper))
            else:
                print("orbitper:  {:.4f} --> {:.4f}".
                      format(old_orbitper, self.orbitper))
            if old_doppzero == -999:
                print("doppzero:  {:.7f} (added)".format(self.doppzero))
            else:
                print("doppzero:  {:.7f} --> {:.7f}".
                      format(old_doppzero, self.doppzero))
            if old_doppmag == -999:
                print("doppmag:   {:.6f} (added)".format(self.doppmag))
            else:
                print("doppmag:   {:.6f} --> {:.6f}".
                      format(old_doppmag, self.doppmag))
            if old_doppmag_v == -999:
                print("doppmagv:  {:.6f} (added)".format(self.doppmag_v))
            else:
                print("doppmagv:  {:.6f} --> {:.6f}".
                      format(old_doppmag_v, self.doppmag_v))

            print("")

if __name__ == "__main__":

    main(sys.argv[1:])
