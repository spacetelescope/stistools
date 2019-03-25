from astropy.io import fits


def initObservation(input, instrument, sci_num):
    """Construct an Observation object for the current mode.

    Parameters
    ----------
    input: str
        The name of an input file.

    instrument: str
        Value of keyword INSTRUME, should be "COS" or "STIS"

    Returns
    -------
    obs: an Observation object
        Information about the observation, mostly from header keywords.
    """

    instrument = instrument.upper()
    if instrument == "STIS":
        obs = Observation(input, sci_num)
    else:
        raise RuntimeError("instrument '{}' is not supported".
                           format(instrument))

    return obs


class Observation(object):
    """Get information about an observation from its headers."""

    def __init__(self, input, sci_ext=1):
        """Invoked by a subclass.

        Parameters
        ----------
        input: str
            The name of an input file.
        """

        self.input = input

        self.sci_ext = sci_ext
        self.ra_targ = None
        self.dec_targ = None
        self.cenwave = None
        self.expstart = None
        self.expend = None
        self.dispersion = None

    def getInfo(self):
        """Get information about the exposure."""

        fd = fits.open(self.input, mode="readonly")

        phdr = fd[0].header
        hdr = fd['sci', self.sci_ext].header

        self.ra_targ = phdr["ra_targ"]
        self.dec_targ = phdr["dec_targ"]
        self.cenwave = phdr.get("cenwave", default=0)

        if self.cenwave <= 0:
            raise ValueError("CENWAVE = %d" % self.cenwave)

        if phdr["detector"] == "CCD":
            highres_factor = 1.0
        else:
            highres_factor = 2.0  # either MAMA detector

        self.expstart = hdr["expstart"]
        self.expend = hdr["expend"]
        cd1_1 = hdr.get("cd1_1", 1.)
        ltm1_1 = hdr.get("ltm1_1", 1.)
        self.dispersion = cd1_1 * ltm1_1 / highres_factor

        if self.dispersion == 0.:
            raise ValueError("dispersion is zero")

        fd.close()
