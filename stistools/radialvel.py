import numpy as N

DEG_RAD = N.pi / 180.                   # degrees to radians
ARCSEC_RAD = N.pi / (180.*3600.)        # arcseconds to radians

REFDATE = 51544.5               # MJD for 2000 Jan 1, 12h UT
KM_AU = 1.4959787e8             # kilometers per astronomical unit
SEC_DAY = 86400.                # seconds per day


def radialVel(ra_targ, dec_targ, mjd):
    """Compute the heliocentric velocity of the Earth.

    This function computes the radial velocity of a target based on the
    Earth's orbital velocity around the Sun.  The space motion of the
    target is not taken into account.  That is, the radial velocity is
    just the negative of the component of the Earth's orbital velocity
    in the direction toward the target.

    Parameters
    -----------
    ra_targ : float
        right ascension of the target (degrees)
    dec_targ : float
        declination of the target (degrees)
    mjd : float
        Modified Julian Date at the time of observation

    Returns
    --------
    radial_vel : float
        the radial velocity in km/s
    """

    # Convert target position to rectangular coordinate unit vector.
    ra = ra_targ * DEG_RAD
    dec = dec_targ * DEG_RAD
    target = N.zeros(3, dtype=N.float64)
    target[0] = N.cos(dec) * N.cos(ra)
    target[1] = N.cos(dec) * N.sin(ra)
    target[2] = N.sin(dec)

    # Precess the target coordinates from J2000 to the current date.
    target = precess(mjd, target)

    # Get the Earth's velocity vector (km/sec).
    velocity = earthVel(mjd)

    # Dot product.
    vel_r = N.dot(velocity, target)

    return -vel_r


def earthVel(mjd):
    """Compute and return the velocity of the Earth at the specified time.

    This function computes the Earth's orbital velocity around the Sun
    in celestial rectangular coordinates.  The expressions are from the
    Astronomical Almanac, p C24, which gives low precision formulas for
    the Sun's coordinates.  We'll apply these formulas directly to get
    the velocity of the Sun relative to the Earth, then we'll convert to
    km per sec and change the sign to get the velocity of the Earth.

    Notes
    -----
    We get the velocity of the Sun relative to the Earth as follows:

    The velocity in the ecliptic plane with the X-axis aligned with the
    radius vector is:

      - Vx = radius_dot,
      - Vy = radius * elong_dot,
      - Vz = 0

    where:

      - radius is the radial distance from Earth to Sun
      - elong is the ecliptic longitude of the Sun
      - eps is the obliquity of the ecliptic
      - _dot means the time derivative

    Rotate in the XY-plane by elong to get the velocity in ecliptic
    coordinates::

      radius_dot * cos (elong) - radius * elong_dot * sin (elong)
      radius_dot * sin (elong) + radius * elong_dot * cos (elong)
      0

    Rotate in the YZ-plane by eps to get the velocity in equatorial
    coordinates::

       radius_dot * cos (elong) - radius * elong_dot * sin (elong)
       (radius_dot * sin (elong) + radius * elong_dot * cos (elong)) * cos (eps)
       (radius_dot * sin (elong) + radius * elong_dot * cos (elong)) * sin (eps)

    Parameters
    ----------
    mjd : float
        time, Modified Julian Date

    Returns
    -------
    vel :  ndarray
        the velocity vector of the Earth around the Sun, in
        celestial coordinates (shape=(3,),ndtype=float64)
    """

    # All angular values are in radians.
    # g      = mean anomaly
    # L      = mean longitude, corrected for aberration
    # elong  = ecliptic longitude
    # radius = distance to Sun in AU
    # eps    = obliquity of ecliptic

    # time in days since JD 2451545.0
    tdays = mjd - REFDATE

    g_dot = 0.9856003 * DEG_RAD
    L_dot = 0.9856474 * DEG_RAD

    eps = (23.439 - 0.0000004 * tdays) * DEG_RAD

    g = ((357.528 + 0.9856003 * tdays) % 360.) * DEG_RAD
    L = ((280.461 + 0.9856474 * tdays) % 360.) * DEG_RAD

    #           1.915 degrees          0.02 degree
    elong = L + 0.033423 * N.sin(g) + 0.000349 * N.sin(2.*g)
    elong_dot = L_dot + \
                0.033423 * N.cos(g) * g_dot + \
                                       0.000349 * N.cos(2.*g) * 2.*g_dot

    radius = 1.00014 - 0.01671 * N.cos(g) - 0.00014 * N.cos(2.*g)
    radius_dot =       0.01671 * N.sin(g) * g_dot + \
                                             0.00014 * N.sin(2.*g) * 2.*g_dot

    x_dot = radius_dot * N.cos(elong) - \
                radius * N.sin(elong) * elong_dot

    y_dot = radius_dot * N.cos(eps) * N.sin(elong) + \
                radius * N.cos(eps) * N.cos(elong) * elong_dot

    z_dot = radius_dot * N.sin(eps) * N.sin(elong) + \
                radius * N.sin(eps) * N.cos(elong) * elong_dot

    # Convert to km/sec with Sun as origin.
    velocity = N.zeros(3, dtype=N.float64)
    velocity[0] = -x_dot * KM_AU / SEC_DAY
    velocity[1] = -y_dot * KM_AU / SEC_DAY
    velocity[2] = -z_dot * KM_AU / SEC_DAY

    return velocity


def precess(mjd, target):
    """Precess target coordinates from J2000 to the date mjd.

    Notes
    -----
    target can be a single vector, e.g. [x0, y0, z0], or it can be
    a 2-D array; in the latter case, the shape should be (n,3)::

        target = [[x0, x1, x2, x3, x4],
                  [y0, y1, y2, y3, y4],
                  [z0, z1, z2, z3, z4]]

    The algorithm used in this function was based on [1]_ and [2]_.

    References
    ----------
    .. [1] Lieske, et al. 1976, Astron & Astrophys vol 58, p 1.

    .. [2] J.H. Lieske, 1979, Astron & Astrophys vol 73, 282-284.

    Parameters
    -----------
    mjd : float
        time, Modified Julian Date
    target : array_like object
        unit vector pointing toward the target, J2000 coordinates

    Returns
    -------
    vector : ndarray
        the target vector (or matrix) precessed to mjd as an
        array object of type float64 and the same shape as target,
        i.e. either (3,) or (n,3)
    """

    target_j2000 = N.array(target, dtype=N.float64)
    target_mjd = target_j2000.copy()

    dt = (mjd - REFDATE) / 36525.
    dt2 = dt**2
    dt3 = dt**3

    zeta = (2306.2181 * dt + 0.30188 * dt2 + 0.017998 * dt3) * ARCSEC_RAD

    z = (2306.2181 * dt + 1.09468 * dt2 + 0.018203 * dt3) * ARCSEC_RAD

    theta = (2004.3109 * dt - 0.42665 * dt2 - 0.041833 * dt3) * ARCSEC_RAD

    cos_zeta = N.cos(zeta)
    sin_zeta = N.sin(zeta)
    cos_z = N.cos(z)
    sin_z = N.sin(z)
    cos_theta = N.cos(theta)
    sin_theta = N.sin(theta)

    # Create the rotation matrix.
    a = N.identity(3, dtype=N.float64)

    a[0, 0] =  cos_z * cos_theta * cos_zeta - sin_z * sin_zeta
    a[0, 1] = -cos_z * cos_theta * sin_zeta - sin_z * cos_zeta
    a[0, 2] = -cos_z * sin_theta

    a[1, 0] =  sin_z * cos_theta * cos_zeta + cos_z * sin_zeta
    a[1, 1] = -sin_z * cos_theta * sin_zeta + cos_z * cos_zeta
    a[1, 2] = -sin_z * sin_theta

    a[2, 0] =  sin_theta * cos_zeta
    a[2, 1] = -sin_theta * sin_zeta
    a[2, 2] =  cos_theta

    # Convert to matrix objects.
    m_a = N.matrix(a)
    m_target_j2000 = N.matrix(target_j2000)

    # The prefix "m_" indicates that the product is actually a matrix.
    m_target_mjd = m_a * m_target_j2000.T

    # Return a simple array (rather than a matrix).
    return m_target_mjd.T.A[0]
