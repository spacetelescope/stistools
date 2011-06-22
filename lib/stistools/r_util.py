from __future__ import division         # confidence high
import os
import os.path
import copy

def expandFileName (filename):
    """Expand environment variable in a file name.

    If the input file name begins with either a Unix-style or IRAF-style
    environment variable (e.g. $lref/name_dqi.fits or lref$name_dqi.fits
    respectively), this routine expands the variable and returns a complete
    path name for the file.

    @param filename:  a file name, possibly including an environment variable
    @type filename:  string

    @return:  the file name with environment variable expanded
    @rtype:  string
    """

    n = filename.find ("$")
    if n == 0:
        if filename != NOT_APPLICABLE:
            # Unix-style file name.
            filename = os.path.expandvars (filename)
    elif n > 0:
        # IRAF-style file name.
        temp = "$" + filename[0:n] + os.sep + filename[n+1:]
        filename = os.path.expandvars (temp)
        # If filename contains "//", delete all of them.
        double_sep = os.sep + os.sep
        filename = filename.replace(double_sep,os.sep)

    return filename

def interpolate (x, values, xp):
    """Interpolate.

    Linear interpolation is used.  If the specified indpendent variable
    value xp is outside the range of the array x, the first (or last)
    value in values will be returned.

    @param x:  array of independent variable values
    @type x:  a sequence object, e.g. an array, int or float

    @param values:  array of dependent variable values
    @type values:  a sequence object, e.g. an array (not character)

    @param xp:  independent variable value at which to interpolate
    @type xp:  int or float

    @return:  linearly interpolated value
    @rtype:  the same type as one element of values
    """

    nvalues = len (values)

    if nvalues == 1 or xp <= x.item(0):
        value = copy.deepcopy (values[0])
    elif xp >= x.item(nvalues-1):
        value = copy.deepcopy (values[nvalues-1])
    else:
        # search for independent variable values that bracket the specified xp
        for i in range (nvalues-1):
            x0 = x.item(i)
            x1 = x.item(i+1)
            if xp >= x0 and xp <= x1:
                if x0 == x1:
                    value = copy.deepcopy (values[i])
                else:
                    p = float (x1 - xp)
                    q = float (xp - x0)
                    value = (p * values[i] + q * values[i+1]) / (p + q)
                break

    return value
