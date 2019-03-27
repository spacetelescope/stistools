#!/usr/bin/env python
import numpy as np
from astropy.io import fits
import os.path
from scipy import signal
from scipy import ndimage as ni

from stsci.tools import gfit, linefit
from stsci.tools import fileutil as fu

__doc__ = """
Refine a STIS trace table.

- A trace is generated from the science file and a trace
  center is computed.
- The two traces bracketing the trace center are extracted
  from the trace table and interpolated
- The correction is computed as the difference between the
  linear fit to the science and interpolated traces
- The correction is applied to all traces in the trace file
  for that particular OPT_ELEM and CENWAVE
- A new trace table is written to the current directory and
  the relevant keywords are updates in the header of the input file.

Examples
--------

Simple example of running mktrace on a STIS file named 'file.fits':

>>> import mktrace
>>> mktrace.mktrace('file.fits', [tracecen=509.4], [weights=[(x1,x2),(x3,x4)])


:Authors:

- Author (IDL): Linda Dressel
- Python version: Nadia Dencheva

"""

__version__ = '2.0.0'
__vdate__ = '2017-03-20'


def mktrace(fname, tracecen=0.0, weights=None):
    """
    Refine a stis spectroscopic trace.
    """
    try:
        hdulist = fits.open(fname)
    except IOError:
        print("\nUNABLE TO OPEN FITS FILE: {} \n".format(fname))
        return

    data = hdulist[1].data
    hdr0 = hdulist[0].header
    hdr1 = hdulist[1].header
    hdulist.close()

    kwinfo = getKWInfo(hdr0, hdr1)
    if kwinfo['instrument'] != 'STIS':
        print("This trace tool works only on STIS spectroscopic observations.\n")
        print("Not processing file {}.\n".format(fname))
        return

    sizex, sizey = data.shape
    if weights is None:
        wei = np.ones(sizey)
    else:
        if not iterable(weights) or not iterable(weights[0]):
            print("Weights must be a list of tuples, for example:\n")
            print("weights=[(23, 45),(300, 670)] \n")
            return
        wei = np.zeros(sizey)

        for i in np.arange(len(weights)):
            for j in np.arange(weights[i][0], weights[i][1]):
                wei[j] = 1

    # wind are weights indices in the image frame which may be a subarray
    wind = np.nonzero(wei)[0]

    tr = Trace(fname, kwinfo)
    a2center, trace1024 = tr.generateTrace(data, kwinfo, tracecen=tracecen, wind=wind)
    # compute the full frame a2center
    ffa2center = a2center*kwinfo['binaxis2']
    tr_ind, a2disp_ind = tr.getTraceInd(ffa2center)
    tr2 = tr.readTrace(tr_ind)
    if tr_ind != a2disp_ind[0]:
        tr1 = tr.readTrace(tr_ind - 1)
        interp_trace = trace_interp(tr1, tr2, ffa2center)
    else:
        interp_trace = tr2

    # convert the weights array into full frame
    ind = np.nonzero(wei)[0] * kwinfo['binaxis1']
    w = np.zeros(1024)
    w[ind] = 1

    X = np.arange(1024).astype(np.float)
    sparams = linefit.linefit(X, trace1024, weights=w)
    rparams = linefit.linefit(X, interp_trace, weights=w)
    sciline = sparams[0] + sparams[1] * X
    refline = rparams[0] + rparams[1] * X

    deltaline = sciline - refline

    # create a complete trace similar to a row in a _1dt file
    # used only for debugging
    tr._a2displ = trace1024
    tr._a1center = tr1['a1center']
    tr._a2center = a2center
    tr._nelem = tr1['nelem']
    tr._pedigree = tr1['pedigree']
    tr._snr_thresh = tr1['snr_thresh']

    tr.writeTrace(fname, sciline, refline, interp_trace,
                  trace1024, tr_ind, a2disp_ind)

    #print 'time', time.time()-start
    #the minus sign is for consistency withthe way x2d reports the rotation
    print("Traces were rotated by {:.10f} degrees \n".format((-(sparams[1]-rparams[1])*180 / np.pi)))
    print('trace is centered on row {:.10f}'.format(tr._a2center))
    return tr


def iterable(v):
    try:
        len(v)
        return True
    except TypeError:
        return False


def interp(y, n):
    """
    Given a 1D array of size m, interpolates it to a size n (m < n).
    """
    m = float(len(y))
    x = np.arange(m)
    i = np.arange(n,dtype=np.float)
    xx = i * (m-1)/n
    xind = np.searchsorted(x, xx)-1
    yy = y[xind]+(xx-x[xind])*(y[xind+1]-y[xind])/(x[xind+1]-x[xind])

    return yy


def trace_interp(tr1, tr2, cen):

    a2disp1 = tr1['a2displ']
    a2disp2 = tr2['a2displ']
    za2disp1 = a2disp1 - a2disp1[512]
    za2disp2 = a2disp2 - a2disp2[512]
    high = tr2['a2center']
    low = tr1['a2center']
    n2 = (cen - low) / (high - low)
    n1 = 1.0 - n2
    interp_trace = n1 * za2disp1 + n2 * za2disp2

    return interp_trace


def getKWInfo(hdr0, hdr1):
    kwinfo = {}
    kwinfo['instrument'] = hdr0['INSTRUME']
    kwinfo['detector'] = hdr0['DETECTOR']
    if kwinfo['detector'] == "CCD":
        kwinfo['binaxis2'] = hdr0['BINAXIS2']
        kwinfo['binaxis1'] = hdr0['BINAXIS1']
    else:
        kwinfo['binaxis2'] = 1
        kwinfo['binaxis1'] = 1
    kwinfo['crpix2'] = hdr1['CRPIX2']
    kwinfo['ltv2'] = hdr1['LTV2']
    kwinfo['sizaxis2'] = hdr0['sizaxis2']
    kwinfo['opt_elem'] = hdr0['OPT_ELEM']
    kwinfo['cenwave'] = hdr0['CENWAVE']
    kwinfo['sporder'] = hdr1['SPORDER']
    kwinfo['sptrctab'] = hdr0['SPTRCTAB']

    return kwinfo


class Trace:
    """ Trace class for a crj or flt file.

    Notes
    -----
    tr=Trace(file)
    file is a crj or flt file.

    opt_elem, cenwave, sporder are read from the header of the science file
    a2center is a2center of the trace generated from the science file

    tr_ind= tr.getTraceInd(a2center)

    tr_ind is the index of the row in the trace file which brackets
    from below a2center as computed fro the generated trace

    tr.readTrace(tr_ind)

    a2center = tr.generateTrace(...)

    """
    def __init__(self, file, kwinfo):
        self._opt_elem = kwinfo['opt_elem']
        self._cenwave = kwinfo['cenwave']
        self._sporder = kwinfo['sporder']
        self._nelem = None
        self._a2displ = None
        self._a1center = None
        self._a2center = None
        self._snr_thresh = None
        self._pedigree = None
        self.sptrctabname = kwinfo['sptrctab']
        self.sptrctab = self.openTraceFile(fu.osfn(self.sptrctabname))

    def openTraceFile(self, filename):

        """
        Returns a spectrum trace table
        """
        if filename is not None:
            try:
                f = fits.open(filename)
            except IOError:
                print("Could not open file {}.\n".format(filename))
                return

            tab = f[1].data
            f.close()
            return tab
        else:
            print("A valid 1-D SPECTRUM TRACE TABLE is required.\n")
            return None

    def getTraceInd(self, a2center):
        """
        Finds the first trace in the trace table whose A2CENTER is larger
        than the specified a2center.
        """
        opt_ind = self.sptrctab.field('OPT_ELEM') == self._opt_elem
        cen_ind = self.sptrctab.field('CENWAVE') == self._cenwave
        sp_ind = self.sptrctab.field('SPORDER') == self._sporder
        a2disp_ind = opt_ind & cen_ind & sp_ind
        ind = np.nonzero(a2disp_ind)
        i = np.nonzero(self.sptrctab[ind].field('A2CENTER') > a2center)[0][0] + ind[0][0]

        return i, a2disp_ind

    def readTrace(self, tr_ind):
        """
        reads the specified row from the 1dttab.fits
        """

        tr = {}
        tr['nelem'] = self.sptrctab[tr_ind].field('NELEM')
        tr['a2displ'] = self.sptrctab[tr_ind].field('A2DISPL')
        tr['a1center'] = self.sptrctab[tr_ind].field('A1CENTER')
        tr['a2center'] = self.sptrctab[tr_ind].field('A2CENTER')
        tr['snr_thresh'] = self.sptrctab[tr_ind].field('SNR_THRESH')
        tr['pedigree'] = self.sptrctab[tr_ind].field('PEDIGREE')

        return tr

    def writeTrace(self, fname, sciline, refline, interp_trace, trace1024,
                   tr_ind, a2disp_ind):

        """
        The 'writeTrace' method performs the following steps:

          - Adds sciline-refline to all traces with the relevent OPT_ELEM,
            CENWAVE and SPORDER.
          - Writes the new trace table to the current directory.
          - Updates the SPTRCTAB keyword in the header to point to the new table.
          - Writes out fits files with the

            - science trace - '_sci'
            - the fit to the science trace - '_scifit'
            - the interpolated trace - '_interp'
            - the linear fit to the interpolated trace - '_interpfit'

        """
        fpath = fu.osfn(self.sptrctabname)
        infile = fname.split('.')
        newname = infile[0] + '_1dt.' + infile[1]

        # refine all traces for this CENWAVE, OPT_ELEM
        fu.copyFile(fpath, newname)
        hdulist = fits.open(newname, mode='update')
        tab = hdulist[1].data
        ind = np.nonzero(a2disp_ind)[0]
        for i in np.arange(ind[0], ind[-1] + 1):
            tab[i].setfield('A2DISPL', tab[i].field('A2DISPL') + (sciline - refline))
        if 'DEGPERYR' in tab.names:
            for i in np.arange(ind[0], ind[-1] + 1):
                tab[i].setfield('DEGPERYR', 0.0)

        hdulist.flush()
        hdulist.close()

        # update SPTRCTAB keyword in the science file primary header
        hdulist = fits.open(fname, mode='update')
        hdr0 = hdulist[0].header
        hdr0['SPTRCTAB'] = newname
        hdulist.close()

        # write out the fit to the interpolated trace ('_interpfit' file)
        refhdu = fits.PrimaryHDU(refline)
        refname = infile[0] + '_1dt_interpfit.' + infile[1]
        if os.path.exists(refname):
            os.remove(refname)
        refhdu.writeto(refname)

        # write out the interpolated trace ('_interp' file)
        inthdu = fits.PrimaryHDU(interp_trace)
        intname = infile[0] + '_1dt_interp.' + infile[1]
        if os.path.exists(intname):
            os.remove(intname)
        inthdu.writeto(intname)

        # write out the the fit to the science trace ('_scifit' file)
        scihdu = fits.PrimaryHDU(sciline)
        sciname = infile[0] + '_1dt_scifit.' + infile[1]
        if os.path.exists(sciname):
            os.unlink(sciname)
        scihdu.writeto(sciname)

        # write out the science trace ('_sci' file)
        trhdu = fits.PrimaryHDU(trace1024)
        trname = infile[0] + '_1dt_sci.' + infile[1]
        if os.path.exists(trname):
            os.unlink(trname)
        trhdu.writeto(trname)

    def generateTrace(self, data, kwinfo, tracecen=0.0, wind=None):
        """
        Generates a trace from a science file.
        """
        if kwinfo['sizaxis2'] is not None and kwinfo['sizaxis2'] < 1023:
            subarray = True
        else:
            subarray = False

        if tracecen == 0:
            if subarray:
                _tracecen = kwinfo['sizaxis2'] / 2.0
            else:
                _tracecen = kwinfo['crpix2']
        else:
            _tracecen = tracecen

        sizex, sizey = data.shape
        subim_size = 40
        y1 = int(_tracecen - subim_size/2.)
        y2 = int(_tracecen + subim_size/2.)
        if y1 < 0:
            y1 = 0
        if y2 > (sizex -1):
            y2 = sizex - 1
        specimage = data[y1:y2+1, :]
        smoytrace = self.gFitTrace(specimage, y1, y2)
        yshift = int(np.median(smoytrace) - 20)
        y1 = y1 + yshift
        y2 = y2 + yshift
        if y1 < 0:
            y1 = 0
        if y2 > sizex:
            y2 = sizex
        specimage = data[y1:y2+1, :]
        smoytrace = self.gFitTrace(specimage, y1, y2)
        med11smoytrace = ni.median_filter(smoytrace, 11)
        med11smoytrace[0] = med11smoytrace[2]
        diffmed = abs(smoytrace - med11smoytrace)
        tolerence = 3 * np.median(abs(smoytrace[wind] - med11smoytrace[wind]))
        if tolerence < 0.1:
            tolerence = 0.1
        badpoint = np.where(diffmed > tolerence)[0]
        if len(badpoint) != 0:
            np.put(smoytrace, badpoint, med11smoytrace[badpoint])

        # Convolve with a gaussian to smooth it.
        fwhm = 10.
        sigma = fwhm / 2.355
        gaussconvxsmoytrace = ni.gaussian_filter1d(smoytrace, sigma)

        # Compute the trace center as the median of the pixels
        # with nonzero weights.
        tracecen = np.median(gaussconvxsmoytrace[wind])
        gaussconvxsmoytrace = gaussconvxsmoytrace - tracecen
        trace1024 = interp(gaussconvxsmoytrace, 1024) * kwinfo['binaxis2']
        tracecen = tracecen + y1 + 1.0
        if subarray:
            tracecen = tracecen - kwinfo['ltv2']
        self.trace1024 = trace1024
        return tracecen, trace1024

    def gFitTrace(self, specimage, y1, y2):
        """
        Fit a gaussian to each column of an image.
        """

        sizex, sizey = specimage.shape
        smoytrace = np.zeros(sizey).astype(np.float)
        boxcar_kernel = signal.boxcar(3) / 3.0

        for c in np.arange(sizey):
            col = specimage[:, c]
            col = col - np.median(col)
            smcol = ni.convolve(col, boxcar_kernel).astype(np.float)
            fit = gfit.gfit1d(smcol, quiet=1, maxiter=15)
            smoytrace[c] = fit.params[1]

        return np.array(smoytrace)
