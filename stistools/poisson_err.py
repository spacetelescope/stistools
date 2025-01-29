#! /usr/bin/env python
import numpy as np
from astropy.stats import poisson_conf_interval
from astropy.io import fits

__doc__ = """
The task :func:`poisson_err` calculates errors according to a Poisson
confidence interval using astropy.stats.poisson_conf_interval. By default, the
pipeline assumes that the errors follow a root-N Gaussian appoximation. This
re-calculation can be important for datasets that have low detector counts
(e.g., FUV continuum).

The input file for :func:`poisson_err` is a STIS x1d file. :func:`poisson_err`
reads in the x1d file, calculates the correct 1-sigma upper and lower Poisson
confidence interval limits, and creates new columns in the fits file.

Because we use NET_ERROR to back out the pre-dark-subtracted counts, we are
assuming that err=sqrt(I), which is valid only for MAMA observations. CCD
observations will also have gain, bias, and readnoise corrections.

File will be written to the filename provided by "output", creating a new file
or overwriting an existing one.

Examples
--------

:func:`poisson_err` with default values:

>>> import stistools
>>> stistools.poisson_err.poisson_err("oep502040_x1d.fits",
...                                 "oep502040_poisson_x1d.fits")
"""

__taskname__ = "poisson_err"
__version__ = "1.0"
__vdate__ = "23-August-2024"
__author__ = "J. Lothringer"


def poisson_err(x1dfile, output,  verbose=True):
    """Adds Poisson Confidence Interval columns to extracted spectra.
       For use with MAMA data only. Works with echelle or first-order
       spectra.

        Parameters
        ----------
        x1dfile: str
            Input MAMA x1d file from which we will calculate the poisson
            confidence interval.

        output: str
            Name of the output FITS file. Will overwrite existing file.

        verbose: bool
            Prints a completion message before function returns.

        Returns
        -------
        Nothing is returned directly, but the file is written to output.
    """
    # Open the file
    with fits.open(x1dfile) as hdu:
        # Warning if using CCD data
        try:
            instrument = hdu[0].header['INSTRUME'].strip()
            detector = hdu[0].header['DETECTOR'].strip()
            if (instrument != 'STIS') or \
               (detector not in {'NUV-MAMA', 'FUV-MAMA'}):
                raise ValueError
        except (ValueError, IndexError,) as e:
            raise ValueError(f"{x1dfile}: Not a STIS/MAMA observation.") from e

        for i in range(1, len(hdu)):
            # Read in important rows and numbers
            exptime = hdu[i].header['EXPTIME']
            flux = hdu[i].data['FLUX']
            net = hdu[i].data['NET']
            gross = hdu[i].data['GROSS']
            net_err = hdu[i].data['NET_ERROR']

            # Set up our arrays to loop through the orders
            lo = np.zeros(gross.shape, dtype=np.float32)
            up = np.zeros(gross.shape, dtype=np.float32)
            Ns = np.zeros(gross.shape, dtype=np.float32)
            for order in range(gross.shape[0]):
                # Rather than recalculate N from gross counts, let us
                # use NET_ERROR to go backward to N.
                # That way, dark counts and everything are preserved!
                # But we do still have to convert from C/s to counts.
                N = (net_err[order]*exptime)**2

                # Set negative counts to zero and save
                N[N < 0.0] = 0.0
                Ns[order] = N

                # Calculate PCI (Poisson Confidence Interval)
                lo[order], up[order] = \
                    poisson_conf_interval(N, interval='frequentist-confidence')
                # poisson_conf_interval gives us the interval,
                # not the 1-sigma size, so calculate it
                lo[order] = N - lo[order]
                up[order] = up[order] - N
            # Convert back to count rate
            lo_rate = lo / exptime
            up_rate = up / exptime

            # And convert to flux
            lo_flux = (lo / exptime) * (flux / net)
            up_flux = (up / exptime) * (flux / net)

            # Create new columns, matching format and units
            dim_str = str(gross.shape[1])
            col_lo = fits.Column(name='NET_ERROR_PCI_LOW',
                                 array=lo_rate, format=dim_str+'E',
                                 unit='Counts/s', disp='G15.7')

            col_up = fits.Column(name='NET_ERROR_PCI_UP',
                                 array=up_rate, format=dim_str+'E',
                                 unit='Counts/s', disp='G15.7')

            col_flux_lo = fits.Column(name='ERROR_PCI_LOW',
                                      array=lo_flux, format=dim_str+'E',
                                      unit='erg/s/cm**2/Angstrom',
                                      disp='G15.7')

            col_flux_up = fits.Column(name='ERROR_PCI_UP',
                                      array=up_flux, format=dim_str+'E',
                                      unit='erg/s/cm**2/Angstrom',
                                      disp='G15.7')

            # Adding in the Ns in case one wants to verify or use
            # a different Poisson interval approximation
            col_Ns = fits.Column(name='N_COUNTS', array=Ns,
                                 format=dim_str+'E', unit='Counts',
                                 disp='G15.7')

            # Append new columns to the existing data
            new_cols = fits.ColDefs([col_lo, col_up,
                                     col_flux_lo, col_flux_up, col_Ns])
            hdu[i] = fits.BinTableHDU.from_columns(
                hdu[i].columns + new_cols, header=hdu[i].header)

        # Save the changes
        hdu.writeto(output, overwrite=True)
        if verbose:
            print(
                f"Added ERROR_PCI_LOW, ERROR_PCI_UP, "
                f"NET_ERROR_PCI_LOW, NET_ERROR_PCI_UP, and N_counts "
                f"columns to {x1dfile} in {output}"
            )
        return
