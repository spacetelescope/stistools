import os

import numpy as np
from astropy.io import fits

def initObservation(input, instrument):
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
    if instrument == "COS":
        obs = COSObs(input)
    elif instrument == "STIS":
        obs = STISObs(input)
    else:
        raise RuntimeError("instrument '{}' is not supported".format(instrument))

    return obs

class Observation(object):
    """Get information about an observation from its headers."""

    def __init__(self, input):
        """Invoked by a subclass.

        Parameters
        ----------
        input: str
            The name of an input file.
        """

        self.input = input

        self.ra_targ = None
        self.dec_targ = None
        self.cenwave = None
        self.expstart = None
        self.expend = None
        self.dispersion = None

class COSObs(Observation):
    """A COS observation."""

    def __init__(self, input):
        """Constructor.

        Parameters
        ----------
        input: str
            The name of an input file.
        """

        Observation.__init__(self, input)

    def getInfo(self):
        """Get information about the exposure."""

        fd = fits.open(self.input, mode="readonly")

        phdr = fd[0].header
        hdr = fd[1].header

        self.ra_targ = phdr["ra_targ"]
        self.dec_targ = phdr["dec_targ"]
        self.cenwave = phdr.get("cenwave", default=0)
        exptype = phdr.get("exptype", default="missing")
        if self.cenwave <= 0:
            opt_elem = phdr.get("opt_elem", default="missing")
            targname = phdr.get("targname", default="missing")
            if opt_elem.startswith("MIR"):
                message = "OPT_ELEM = %s" % opt_elem
            elif targname == "DARK":
                message = "This is a DARK exposure."
            elif exptype.startswith("ACQ"):
                message = "EXPTYPE = %s" % exptype
            else:
                message = "CENWAVE = %d" % self.cenwave
            raise ValueError(message)
        if exptype == "WAVECAL":
            print("Warning:  This is a wavecal exposure.")

        self.expstart = hdr["expstart"]
        self.expend = hdr["expend"]

        if hdr["xtension"] == "BINTABLE":
            cd1_1 = hdr.get("tc2_2", 0.)
        elif hdr["xtension"] == "IMAGE":
            cd1_1 = hdr.get("cd1_1", 0.)
        else:
            cd1_1 = 0.
        if cd1_1 == 0. or cd1_1 == 1.:
            # CD matrix is not valid; get dispersion from disptab.
            detector = phdr["detector"]
            aperture = phdr["aperture"]
            if aperture != "BOA":
                aperture = "PSA"
            opt_elem = phdr["opt_elem"]
            fpoffset = phdr.get("fpoffset", "N/A")
            disptab = self.expandFileName(phdr["disptab"])
            if detector == "FUV":
                segment = phdr["segment"]
                middle = 16384. / 2
            else:
                segment = "NUVB"
                middle = 1024. / 2
            filter = {"segment": segment, "opt_elem": opt_elem,
                      "cenwave": self.cenwave, "aperture": aperture}
            if self.findColumn(disptab, "fpoffset"):
                filter["fpoffset"] = fpoffset
            disp_info = self.getTable(disptab, filter)
            if disp_info is None:
                raise RuntimeError("can't get dispersion from disptab")
            ncoeff = disp_info.field("nelem")[0]
            coeff = disp_info.field("coeff")[0][0:ncoeff]
            x = middle
            sum = (ncoeff-1.) * coeff[ncoeff-1]
            for n in range(ncoeff-2, 0, -1):
                sum = sum * x + n * coeff[n]
            self.dispersion = sum
        else:
            self.dispersion = cd1_1

        if self.dispersion == 0.:
            raise ValueError("dispersion is zero")

        fd.close()

    def expandFileName(self, filename):
        """Expand environment variable in a file name.

        If the input file name begins with either a Unix-style or
        IRAF-style environment variable (e.g. $lref/name_dqi.fits or
        lref$name_dqi.fits respectively), this routine expands the
        variable and returns a complete path name for the file.

        Parameters
        ----------
        filename: str
            A name that might include an environment variable.

        Returns
        -------
        filename: str
            The name with the environment variable expanded.
        """

        n = filename.find("$")
        if n == 0:
            if filename != "N/A":
                # Unix-style file name.
                filename = os.path.expandvars(filename)
        elif n > 0:
            # IRAF-style file name.
            temp = "$" + filename[0:n] + os.sep + filename[n+1:]
            filename = os.path.expandvars(temp)
            # If filename contains "//", delete one of them.
            i = filename.find(os.sep + os.sep)
            if i != -1:
                filename = filename[:i+1] + filename[i+2:]

        return filename

    def findColumn(self, table, colname):
        """Return True if colname is found (case-insensitive) in table.

        Parameters
        ----------
        table: str (if name of table) or FITS record object
            Name of table, or data block for a FITS table.

        colname: str
            Name to test for existence in `table`.

        Returns
        -------
        flag: boolean
            True if `colname` is in `table` (without regard to case).
        """

        if type(table) is str:
            fd = fits.open(table, mode="readonly")
            fits_rec = fd[1].data
            fd.close()
        else:
            fits_rec = table

        names = []
        for name in fits_rec.names:
            names.append(name.lower())

        if colname.lower() in names:
            return True
        else:
            return False

    def getTable(self, table, filter):
        """Return the row of a table that matches the filter.

        Parameters
        ----------
        table: str
            Name of the reference table.

        filter: dictionary
            Info for selecting the appropriate row to use; each key in
            the filter is a column name, and if the value in that column
            matches the filter value for some row, that row will be
            included in the set that is returned.

        Returns
        -------
        newdata: FITS record array
            The row or rows that match the filter; in this case there
            should just be one row.
        """

        fd = fits.open(table, mode="readonly")
        data = fd[1].data

        # There will be one element of select_arrays for each non-trivial
        # selection criterion.  Each element of select_arrays is an array
        # of flags, true if the row matches the criterion.
        select_arrays = []
        for key in filter.keys():
            column = data.field(key)
            selected = (column == filter[key])
            select_arrays.append(selected)

        if len(select_arrays) > 0:
            selected = select_arrays[0]
            for sel_i in select_arrays[1:]:
                 selected = np.logical_and(selected, sel_i)
            newdata = data[selected]
        else:
            newdata = fd[1].data.copy()

        fd.close()

        nselect = len(newdata)
        if nselect < 1:
            newdata = None

        return newdata

class STISObs(Observation):
    """A STIS observation."""

    def __init__(self, input):
        """Constructor.

        Parameters
        ----------
        input: str
            The name of an input file.
        """

        Observation.__init__(self, input)

    def getInfo(self):
        """Get information about the exposure."""

        fd = fits.open(self.input, mode="readonly")

        phdr = fd[0].header
        hdr = fd[1].header

        self.ra_targ = phdr["ra_targ"]
        self.dec_targ = phdr["dec_targ"]
        self.cenwave = phdr.get("cenwave", default=0)
        if self.cenwave <= 0:
            raise ValueError("CENWAVE = %d" % self.cenwave)
            
        if phdr["detector"] == "CCD":
            highres_factor = 1.0
        else:
            highres_factor = 2.0                # either MAMA detector

        self.expstart = hdr["expstart"]
        self.expend = hdr["expend"]
        cd1_1 = hdr.get("cd1_1", 1.)
        ltm1_1 = hdr.get("ltm1_1", 1.)
        self.dispersion = cd1_1 * ltm1_1 / highres_factor

        if self.dispersion == 0.:
            raise ValueError("dispersion is zero")

        fd.close()
