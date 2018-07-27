import math
import numpy as np
from astropy.io import fits

STRING_WILDCARD = "ANY"
INT_WILDCARD = -1


def getTable(table, filter, sortcol=None,
             exactly_one=False, at_least_one=False):
    """Return row(s) of a table that match the filter.

    Rows that match every item in the filter (a dictionary of
    column_name=value) will be returned.  If the value in the table is
    STRING_WILDCARD or INT_WILDCARD (depending on the data type of the
    column), that value is considered to match the filter for that column.
    Also, for a given filter key, if the corresponding value in the filter
    is STRING_WILDCARD, the test on filter will be skipped for that key
    (i.e. a wildcard filter element matches any row).

    If more than one row matches the filter, there is an option to sort
    these rows based on the values of one of the table columns.

    It is an error if exactly_one or at_least_one is True but no row
    matches the filter.  A warning will be printed if exactly_one is True
    but more than one row matches the filter.

    Parameters
    -----------
    table : string
        name of the reference table
    filter : dict
        each key is a column name, and the corresponding value
        is a possible table value in that column
    sortcol : string
        the name of a column on which to sort the table rows
        (if there is more than one matching row), or None to disable sorting
    exactly_one : bool
        set this to True if there must be one and only one
        matching row
    at_least_one : bool
        set this to True if there must be at least one
        matching row

    Returns
    -------
    match_rows : rec_array
        an array of the rows of the table that match the filter;
        note that if only one row matches the filter, the function value
        will still be an array

    """

    fd = fits.open(table, mode="readonly")
    data = fd[1].data

    # There will be one element of select_arrays for each non-trivial
    # selection criterion.  Each element of select_arrays is an array
    # of flags, true if the row matches the criterion.
    select_arrays = []
    for key in filter.keys():

        if filter[key] == STRING_WILDCARD:
            continue
        column = data.field(key)
        if len(column) == 0:
            return None
        selected = (column == filter[key])

        # Test for for wildcards in the table.
        wild = None
        if isinstance(column, np.chararray):
            wild = (column == STRING_WILDCARD)
        elif isinstance(column[0], np.integer):
            wild = (column == INT_WILDCARD)
        if wild is not None:
            selected = np.logical_or(selected, wild)

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

    if (exactly_one or at_least_one) and nselect < 1:
        message = "Table has no matching row;\n" + \
                  "table name is " + table + "\n" + \
                  "row selection is " + repr(filter)
        raise RuntimeError(message)

    if exactly_one and nselect > 1:
        print("Table has more than one matching row;")
        print("table name is", table)
        print("row selection is", repr(filter))
        print("only the first will be used.")

    if len(newdata) > 1 and sortcol is not None:
        newdata = sortrows(newdata, sortcol)

    return newdata


def sortrows(rowdata, sortcol, ascend=True):
    """Return a copy of rowdata, sorted on sortcol."""

    if len(rowdata) <= 1:
        return rowdata

    column = rowdata.field(sortcol)
    index = column.argsort()
    if not ascend:
        ind = list(index)
        ind.reverse()
        index = np.array(ind)

    return rowdata[index]


def rotateTrace(trace_info, expstart):
    """Rotate a2displ, if MJD and DEGPERYR are in the trace table.

    Parameters
    -----------
    trace_info : rec_array
        an array of the relevant rows of the table;
        the A2DISPL column will be modified in-place if the MJD and
        DEGPERYR columns are present
    expstart : float
        exposure start time (MJD)
    """

    if expstart < 0:
        return

    # If these columns are not in the table, just return.
    names = []
    for name in trace_info.names:
        names.append(name.lower())
    if "degperyr" in names and "mjd" in names:
        degperyr = trace_info.field("degperyr")
        mjd = trace_info.field("mjd")
    else:
        return

    a2displ = trace_info.field("a2displ")
    nelem = trace_info.field("nelem")
    for i in range(len(trace_info)):
        angle = (degperyr[i] * (expstart - mjd[i]) / 365.25)
        tan_angle = math.tan(angle * math.pi / 180.)
        x = np.arange(nelem[i], dtype=np.float64)
        x -= (nelem[i] // 2)
        a2displ[i][:] -= (x * tan_angle)
