import os
import os.path
import copy


NOT_APPLICABLE = 'n/a'

def expandFileName(filename):
    """Expand environment variable in a file name.

    If the input file name begins with either a Unix-style or IRAF-style
    environment variable (e.g. $lref/name_dqi.fits or lref$name_dqi.fits
    respectively), this routine expands the variable and returns a complete
    path name for the file.

    Parameters
    -----------
    filename : str
        A file name, possibly including an environment variable.

    Returns
    --------
    fullname : str
        The file name with environment variable expanded.
    """

    n = filename.find("$")
    if n == 0:
        if filename != NOT_APPLICABLE:
            # Unix-style file name.
            filename = os.path.expandvars(filename)
    elif n > 0:
        # IRAF-style file name.
        temp = "$" + filename[0:n] + os.sep + filename[n+1:]
        filename = os.path.expandvars(temp)
        # If filename contains "//", delete all of them.
        double_sep = os.sep + os.sep
        filename = filename.replace(double_sep, os.sep)

    return filename


def interpolate(x, values, xp):
    """Interpolate.

    Linear interpolation is used.  If the specified indpendent variable
    value xp is outside the range of the array x, the first (or last)
    value in values will be returned.

    Parameters
    -----------
    x : a sequence object, e.g. an array, int or float
        Array of independent variable values.
    values :  a sequence object, e.g. an array (not character)
        Array of dependent variable values.
    xp : int or float
        Independent variable value at which to interpolate.

    Returns
    -------
    interp_vals : the same type as one element of values
        Linearly interpolated value.

    """

    nvalues = len(values)

    if nvalues == 1 or xp <= x.item(0):
        value = copy.deepcopy(values[0])
    elif xp >= x.item(nvalues-1):
        value = copy.deepcopy(values[nvalues-1])
    else:
        # search for independent variable values that bracket the specified xp
        for i in range(nvalues-1):
            x0 = x.item(i)
            x1 = x.item(i+1)
            if xp >= x0 and xp <= x1:
                if x0 == x1:
                    value = copy.deepcopy(values[i])
                else:
                    p = float(x1 - xp)
                    q = float(xp - x0)
                    value = (p * values[i] + q * values[i+1]) / (p + q)
                break

    return value
