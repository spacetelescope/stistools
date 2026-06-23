#! /usr/bin/env python

import warnings
from copy import deepcopy
from functools import reduce
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
from astropy.table import Table
from .calculate_sensitivity import calculate_sensitivity

__doc__ = r"""
The ``splice`` module concatenates the several orders contained in a STIS
Echelle ``_x1d`` spectrum while co-adding the overlapping sections. This code
emulates the splice module previously implemented on the STSDAS IRAF package.

The function :func:`splice` takes as input an ``_x1d`` Echelle spectrum
(including path) and outputs an Astropy Table containing the wavelength, flux,
flux uncertainty and data-quality flags. The spliced spectrum can be saved as
an additional extension to the ``_x1d`` file by setting ``update_fits`` to
``True``. The spliced spectrum can also be exported to an ascii file by setting
a path and filename to ``output_file``.

The codebase from the ULLYSES project, which is also sponsored by STScI, possess
routines that can be used to co-add and splice STIS Echelle spectra. The main
differences between these codes are the following: 1) ``splice`` does not change
the spectra in non-overlap regions; 2) ``splice`` takes into account
data-quality (DQ) flags and assigns specific DQ values to co-added pixels;
3) ``splice`` only works for STIS spectra and does not concatenate or co-add
COS spectra or multiple STIS spectra.

Examples
--------

Read a spectrum with :func:`read_spectrum` and plot all the different orders
using ``matplotlib`` to visualize how an Echelle extracted spectrum looks like:

.. code-block:: python

   import matplotlib.pyplot as plt
   from stistools import splice
   plt.ion()

   spectrum = splice.read_spectrum('oblh01040_x1d.fits')
   for s in spectrum:
       plt.plot(s['wavelength'], s['flux'],
           linewidth=1, alpha=0.6, color=f"C{s['sporder'] % 2}")
   plt.xlabel(r'Wavelength (Å)')
   plt.ylabel(r'Flux Density (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
   plt.show()

.. image:: splice_eg1.png

Splice the Echelle spectrum orders with :func:`splice` and plot it
using ``matplotlib``:

.. code-block:: python

   plt.clf()

   spliced_spectrum = splice.splice('oblh01040_x1d.fits')
   plt.errorbar(spliced_spectrum['WAVELENGTH'], spliced_spectrum['FLUX'],
                yerr=spliced_spectrum['ERROR'], linewidth=1)
   plt.xlabel(r'Wavelength (Å)')
   plt.ylabel(r'Flux Density (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)')
   plt.show()

.. image:: splice_eg2.png
"""

__taskname__ = "splice"
__version__ = "1.0"
__vdate__ = "23-June-2026"
__author__ = "Leonardo Dos Santos"
__all__ = ["nearest_index", "read_spectrum", "find_overlap", "merge_overlap",
           "concatenate_sections", "splice"]


# Useful tool when dealing with binned data
def nearest_index(array, target_value):
    """
    Finds the index of a value in ``array`` that is closest to ``target_value``.

    Parameters
    ----------
    array : ``numpy.array``
        Target array
    target_value : ``float``
        Target value

    Returns
    -------
    index : ``int``
        Index of the value in ``array`` that is closest to ``target_value``.
    """
    return np.abs(array - target_value).argmin()


# Read the Echelle spectrum based on dataset name and a prefix for file location
def read_spectrum(x1d_input, ext=1, truncate_edge_left=5, truncate_edge_right=5):
    """
    This is a fairly straightforward function to read the spectrum from a `x1d`
    FITS file.

    Parameters
    ----------
    x1d_input : ``str``
        Path and name of the ``*_x1d.fits`` file containing the spectrum.

    ext : ``int``, optional
        FITS extension to read from the X1D file.  Most X1D files contain only
        one SCI extension (i.e., ``REPEATOBS=1``).  Default is ``1``.

    truncate_edge_left : ``int``, optional
        Set the number of low-resolution pixels at the left edge of the detector
        where the spectra should be truncated. If ``0``, then no truncation
        is applied. Default is ``5``.

    truncate_edge_right : ``int``, optional
        Set the number of low-resolution pixels at the right edge of the
        detector where the spectra should be truncated. If ``0``, then no
        truncation is applied. Default is ``5``.

    Returns
    -------
    spectrum : ``list``
        List of all the orders contained in the Echelle spectrum and their
        respective fluxes.
    """
    with fits.open(x1d_input) as hdu:
        header = hdu[0].header
        data = hdu[ext].data
    optical_element = header['OPT_ELEM']
    if optical_element[0] != 'E':
        raise TypeError("This is not an Echelle spectrum.")
    wavelength = data['WAVELENGTH']
    flux = data['FLUX']
    uncertainty = data['ERROR']
    data_quality = data['DQ']
    gross_counts = data['GROSS']
    net_counts = data['NET']
    sporders = data['SPORDER']
    n_orders = len(wavelength)

    # Determine how to truncate the left and right edges with a slice:
    if (not isinstance(truncate_edge_left, int)) or (truncate_edge_left < 0):
        raise ValueError('truncate_edge_left should be a positive integer.')
    if (not isinstance(truncate_edge_right, int)) or (truncate_edge_right < 0):
        raise ValueError('truncate_edge_right should be a positive integer.')
    subset = slice(truncate_edge_left, flux.shape[1] - truncate_edge_right)

    # We index the spectral regions as a `dict` in order to avoid confusion with
    # too many numerical indexes. Also, since the orders are in reverse order,
    # we index them in the opposite way
    spectrum = [{'filename': x1d_input,
                 'ext': ext,
                 'sporder': int(sporders[-i - 1]),
                 'wavelength': wavelength[-i - 1][subset],
                 'flux': flux[-i - 1][subset],
                 'uncertainty': uncertainty[-i - 1][subset],
                 'data_quality': data_quality[-i - 1][subset],
                 'gross': gross_counts[-i - 1][subset],
                 'net': net_counts[-i - 1][subset]}
                for i in range(n_orders)]

    return spectrum


# Identify overlaps in the whole spectrum
def find_overlap(spectrum, extent=1024):
    """
    Find and return the overlapping sections of the Echelle spectrum.

    Parameters
    ----------
    spectrum : ``list``
        List of dictionaries containing the orders of the Echelle spectrum.
        It should resemble the output of ``read_spectrum()``.

    extent : ``int``
        Pixel extent of a full order after truncation.  Typically
        1024 - truncate_left_edge - truncate_right_edge.  Default is
        ``1024``.

    Returns
    -------
    unique_sections : ``list``
        List containing the unique sections of the spectrum.

    overlap_pair_sections : ``list``
        List containing the overlapping pairs of the spectrum.

    overlap_trio_sections : ``list``
        List containing the overlapping trios of the spectrum.
    """
    right_edge = extent - 1

    n_orders = len(spectrum)

    # Identify the wavelength borders of each order
    borders = []
    for order in spectrum:
        borders.append([min(order['wavelength']), max(order['wavelength'])])

    unique_sections = []
    first_trio = []
    second_trio = []
    third_trio = []
    first_pair = []
    second_pair = []

    # The following code is hacky and not pretty to look at, but it works. Sorry
    # for the mess!

    # First we deal with the first orders
    order = spectrum[0]
    wl = order['wavelength']
    idx = [nearest_index(wl, borders[1][0]),
           nearest_index(wl, borders[2][0])]

    # There is always a unique section here
    unique_idx = np.arange(0, idx[0], 1)
    unique_sections.append({'filename': order['filename'],
                            'ext': order['ext'],
                            'sporder': order['sporder'],
                            'wavelength': order['wavelength'][unique_idx],
                            'flux': order['flux'][unique_idx],
                            'uncertainty': order['uncertainty'][unique_idx],
                            'data_quality': order['data_quality'][unique_idx],
                            'gross': order['gross'][unique_idx],
                            'net': order['net'][unique_idx]})

    # There is also a pair overlap. But since it's with the next order, it's
    # considered a "second pair"
    overlap_01 = np.arange(idx[0], idx[1], 1)
    second_pair.append({'filename': order['filename'],
                        'ext': order['ext'],
                        'sporder': order['sporder'],
                        'wavelength': order['wavelength'][overlap_01],
                        'flux': order['flux'][overlap_01],
                        'uncertainty': order['uncertainty'][overlap_01],
                        'data_quality': order['data_quality'][overlap_01],
                        'gross': order['gross'][overlap_01],
                        'net': order['net'][overlap_01]})

    if idx[1] < right_edge:
        # There is a trio overlap
        overlap_012 = np.arange(idx[1], right_edge, 1)
        third_trio.append({'filename': order['filename'],
                           'ext': order['ext'],
                           'sporder': order['sporder'],
                           'wavelength': order['wavelength'][overlap_012],
                           'flux': order['flux'][overlap_012],
                           'uncertainty': order['uncertainty'][overlap_012],
                           'data_quality': order['data_quality'][overlap_012],
                           'gross': order['gross'][overlap_012],
                           'net': order['net'][overlap_012]})

    # Now the second order
    order = spectrum[1]
    wl = order['wavelength']
    idx = [nearest_index(wl, borders[2][0]),
           nearest_index(wl, borders[0][1]),
           nearest_index(wl, borders[3][0])]

    # There are two pairs, potentially one or two trios, and potentially a
    # unique section
    if idx[1] > idx[0]:
        # There is a trio overlap
        overlap_01 = np.arange(0, idx[0], 1)
        overlap_012 = np.arange(idx[0], idx[1], 1)
        unique_1 = None

        if idx[2] < right_edge:
            # There is another trio overlap
            overlap_12 = np.arange(idx[1], idx[2], 1)
            overlap_123 = np.arange(idx[2], right_edge, 1)
        else:
            overlap_123 = None
            overlap_12 = np.arange(idx[0], idx[2], 1)

    else:
        # No trio overlap
        overlap_01 = np.arange(0, idx[1], 1)
        overlap_012 = None
        unique_1 = np.arange(idx[1], idx[0], 1)
        overlap_123 = None
        overlap_12 = np.arange(idx[0], idx[2], 1)

    # Add the to the lists
    if unique_1 is not None:
        unique_sections.append(
            {'filename': order['filename'],
             'ext': order['ext'],
             'sporder': order['sporder'],
             'wavelength': order['wavelength'][unique_1],
             'flux': order['flux'][unique_1],
             'uncertainty': order['uncertainty'][unique_1],
             'data_quality': order['data_quality'][unique_1],
             'gross': order['gross'][unique_1],
             'net': order['net'][unique_1]}
        )

    first_pair.append({'filename': order['filename'],
                       'ext': order['ext'],
                       'sporder': order['sporder'],
                       'wavelength': order['wavelength'][overlap_01],
                       'flux': order['flux'][overlap_01],
                       'uncertainty': order['uncertainty'][overlap_01],
                       'data_quality': order['data_quality'][overlap_01],
                       'gross': order['gross'][overlap_01],
                       'net': order['net'][overlap_01]})
    second_pair.append({'filename': order['filename'],
                        'ext': order['ext'],
                        'sporder': order['sporder'],
                        'wavelength': order['wavelength'][overlap_12],
                        'flux': order['flux'][overlap_12],
                        'uncertainty': order['uncertainty'][overlap_12],
                        'data_quality': order['data_quality'][overlap_12],
                        'gross': order['gross'][overlap_12],
                        'net': order['net'][overlap_12]})

    if overlap_012 is not None:
        second_trio.append({'filename': order['filename'],
                            'ext': order['ext'],
                            'sporder': order['sporder'],
                            'wavelength': order['wavelength'][overlap_012],
                            'flux': order['flux'][overlap_012],
                            'uncertainty': order['uncertainty'][overlap_012],
                            'data_quality': order['data_quality'][overlap_012],
                            'gross': order['gross'][overlap_012],
                            'net': order['net'][overlap_012]})

    if overlap_123 is not None:
        third_trio.append({'filename': order['filename'],
                           'ext': order['ext'],
                           'sporder': order['sporder'],
                           'wavelength': order['wavelength'][overlap_123],
                           'flux': order['flux'][overlap_123],
                           'uncertainty': order['uncertainty'][overlap_123],
                           'data_quality': order['data_quality'][overlap_123],
                           'gross': order['gross'][overlap_123],
                           'net': order['net'][overlap_123]})

    # Now we deal with the third to third before last orders in a loop
    for i in range(n_orders - 4):
        order = spectrum[i + 2]
        wl = order['wavelength']
        idx = [nearest_index(wl, borders[i][1]),
               nearest_index(wl, borders[i + 3][0]),
               nearest_index(wl, borders[i + 1][1]),
               nearest_index(wl, borders[i + 4][0])]
        if idx[0] > 0:
            overlap_idx_012 = np.arange(0, idx[0], 1)
        else:
            overlap_idx_012 = None
        if idx[2] < idx[1]:
            overlap_idx_12 = np.arange(idx[0], idx[2], 1)
            overlap_idx_123 = None
            unique_idx_2 = np.arange(idx[2], idx[1], 1)
            overlap_idx_23 = np.arange(idx[1], idx[3], 1)
        elif idx[2] == idx[1]:
            overlap_idx_12 = np.arange(idx[0], idx[2], 1)
            overlap_idx_123 = None
            unique_idx_2 = np.array([idx[2], ])
            overlap_idx_23 = np.arange(idx[1], idx[3], 1)
        else:
            overlap_idx_12 = np.arange(idx[0], idx[1], 1)
            overlap_idx_123 = np.arange(idx[1], idx[2], 1)
            unique_idx_2 = None
            overlap_idx_23 = np.arange(idx[2], idx[3], 1)
        if idx[3] < right_edge:
            overlap_idx_234 = np.arange(idx[3], right_edge, 1)
        else:
            overlap_idx_234 = None

        if len(overlap_idx_12) > 0:
            first_pair.append(
                {'filename': order['filename'],
                 'ext': order['ext'],
                 'sporder': order['sporder'],
                 'wavelength': order['wavelength'][overlap_idx_12],
                 'flux': order['flux'][overlap_idx_12],
                 'uncertainty': order['uncertainty'][overlap_idx_12],
                 'data_quality': order['data_quality'][overlap_idx_12],
                 'gross': order['gross'][overlap_idx_12],
                 'net': order['net'][overlap_idx_12]}
            )

        if len(overlap_idx_23) > 0:
            second_pair.append(
                {'filename': order['filename'],
                 'ext': order['ext'],
                 'sporder': order['sporder'],
                 'wavelength': order['wavelength'][overlap_idx_23],
                 'flux': order['flux'][overlap_idx_23],
                 'uncertainty': order['uncertainty'][overlap_idx_23],
                 'data_quality': order['data_quality'][overlap_idx_23],
                 'gross': order['gross'][overlap_idx_23],
                 'net': order['net'][overlap_idx_23]}
            )

        if overlap_idx_012 is not None:
            first_trio.append(
                {'filename': order['filename'],
                 'ext': order['ext'],
                 'sporder': order['sporder'],
                 'wavelength': order['wavelength'][overlap_idx_012],
                 'flux': order['flux'][overlap_idx_012],
                 'uncertainty': order['uncertainty'][overlap_idx_012],
                 'data_quality': order['data_quality'][overlap_idx_012],
                 'gross': order['gross'][overlap_idx_012],
                 'net': order['net'][overlap_idx_012]}
            )

        if overlap_idx_123 is not None:
            second_trio.append(
                {'filename': order['filename'],
                 'ext': order['ext'],
                 'sporder': order['sporder'],
                 'wavelength': order['wavelength'][overlap_idx_123],
                 'flux': order['flux'][overlap_idx_123],
                 'uncertainty': order['uncertainty'][overlap_idx_123],
                 'data_quality': order['data_quality'][overlap_idx_123],
                 'gross': order['gross'][overlap_idx_123],
                 'net': order['net'][overlap_idx_123]}
            )

        if overlap_idx_234 is not None:
            third_trio.append(
                {'filename': order['filename'],
                 'ext': order['ext'],
                 'sporder': order['sporder'],
                 'wavelength': order['wavelength'][overlap_idx_234],
                 'flux': order['flux'][overlap_idx_234],
                 'uncertainty': order['uncertainty'][overlap_idx_234],
                 'data_quality': order['data_quality'][overlap_idx_234],
                 'gross': order['gross'][overlap_idx_234],
                 'net': order['net'][overlap_idx_234]}
            )

        if unique_idx_2 is not None:
            unique_sections.append(
                {'filename': order['filename'],
                 'ext': order['ext'],
                 'sporder': order['sporder'],
                 'wavelength': order['wavelength'][unique_idx_2],
                 'flux': order['flux'][unique_idx_2],
                 'uncertainty': order['uncertainty'][unique_idx_2],
                 'data_quality': order['data_quality'][unique_idx_2],
                 'gross': order['gross'][unique_idx_2],
                 'net': order['net'][unique_idx_2]}
            )

    # Now we deal with the last orders. Almost there!
    order = spectrum[-2]
    wl = order['wavelength']
    idx = [nearest_index(wl, borders[-4][1]),
           nearest_index(wl, borders[-1][0]),
           nearest_index(wl, borders[-3][1])]

    # There are two pairs, potentially one or two trios, and potentially a
    # unique section
    if idx[0] > 0:
        # There is a trio overlap
        overlap_012 = np.arange(0, idx[0], 1)
        if idx[2] > idx[1]:
            # There is another trio overlap
            overlap_123 = np.arange(idx[1], idx[2], 1)
            overlap_12 = np.arange(idx[0], idx[1], 1)
            overlap_23 = np.arange(idx[2], right_edge, 1)
            unique_2 = None
        else:
            overlap_123 = None
            unique_2 = np.arange(idx[2], idx[1], 1)
            overlap_12 = np.arange(idx[0], idx[2], 1)
            overlap_23 = np.arange(idx[1], right_edge, 1)
    else:
        overlap_012 = None
        overlap_123 = None
        unique_2 = np.arange(idx[2], idx[1], 1)
        overlap_12 = np.arange(idx[0], idx[2], 1)
        overlap_23 = np.arange(idx[1], right_edge, 1)

    # Add the to the lists
    if unique_2 is not None:
        unique_sections.append(
            {'filename': order['filename'],
             'ext': order['ext'],
             'sporder': order['sporder'],
             'wavelength': order['wavelength'][unique_2],
             'flux': order['flux'][unique_2],
             'uncertainty': order['uncertainty'][unique_2],
             'data_quality': order['data_quality'][unique_2],
             'gross': order['gross'][unique_2],
             'net': order['net'][unique_2]}
        )

    if len(overlap_12) > 0:
        first_pair.append({'filename': order['filename'],
                           'ext': order['ext'],
                           'sporder': order['sporder'],
                           'wavelength': order['wavelength'][overlap_12],
                           'flux': order['flux'][overlap_12],
                           'uncertainty': order['uncertainty'][overlap_12],
                           'data_quality': order['data_quality'][overlap_12],
                           'gross': order['gross'][overlap_12],
                           'net': order['net'][overlap_12]})

    if len(overlap_23) > 0:
        second_pair.append({'filename': order['filename'],
                            'ext': order['ext'],
                            'sporder': order['sporder'],
                            'wavelength': order['wavelength'][overlap_23],
                            'flux': order['flux'][overlap_23],
                            'uncertainty': order['uncertainty'][overlap_23],
                            'data_quality': order['data_quality'][overlap_23],
                            'gross': order['gross'][overlap_23],
                            'net': order['net'][overlap_23]})

    if overlap_012 is not None:
        first_trio.append({'filename': order['filename'],
                           'ext': order['ext'],
                           'sporder': order['sporder'],
                           'wavelength': order['wavelength'][overlap_012],
                           'flux': order['flux'][overlap_012],
                           'uncertainty': order['uncertainty'][overlap_012],
                           'data_quality': order['data_quality'][overlap_012],
                           'gross': order['gross'][overlap_012],
                           'net': order['net'][overlap_012]})

    if overlap_123 is not None:
        second_trio.append({'filename': order['filename'],
                            'ext': order['ext'],
                            'sporder': order['sporder'],
                            'wavelength': order['wavelength'][overlap_123],
                            'flux': order['flux'][overlap_123],
                            'uncertainty': order['uncertainty'][overlap_123],
                            'data_quality': order['data_quality'][overlap_123],
                            'gross': order['gross'][overlap_123],
                            'net': order['net'][overlap_123]})

    # Finally deal with the last order
    order = spectrum[-1]
    wl = order['wavelength']
    idx = [nearest_index(wl, borders[-3][1]),
           nearest_index(wl, borders[-2][1])]

    # There is always a unique section here
    unique_idx = np.arange(idx[1], right_edge, 1)
    unique_sections.append(
        {'filename': order['filename'],
         'ext': order['ext'],
         'sporder': order['sporder'],
         'wavelength': order['wavelength'][unique_idx],
         'flux': order['flux'][unique_idx],
         'uncertainty': order['uncertainty'][unique_idx],
         'data_quality': order['data_quality'][unique_idx],
         'gross': order['gross'][unique_idx],
         'net': order['net'][unique_idx]}
    )

    # There is also a pair overlap. But since it's with the previous order, it's
    # considered a "first pair"
    overlap_23 = np.arange(idx[0], idx[1], 1)

    if len(overlap_23) > 0:
        first_pair.append(
            {'filename': order['filename'],
             'ext': order['ext'],
             'sporder': order['sporder'],
             'wavelength': order['wavelength'][overlap_23],
             'flux': order['flux'][overlap_23],
             'uncertainty': order['uncertainty'][overlap_23],
             'data_quality': order['data_quality'][overlap_23],
             'gross': order['gross'][overlap_23],
             'net': order['net'][overlap_23]}
        )

    if idx[0] > 0:
        # There is a trio overlap
        overlap_123 = np.arange(0, idx[0], 1)
        first_trio.append(
            {'filename': order['filename'],
             'ext': order['ext'],
             'sporder': order['sporder'],
             'wavelength': order['wavelength'][overlap_123],
             'flux': order['flux'][overlap_123],
             'uncertainty': order['uncertainty'][overlap_123],
             'data_quality': order['data_quality'][overlap_123],
             'gross': order['gross'][overlap_123],
             'net': order['net'][overlap_123]}
        )

    # With all that done, we assemble the overlap sections into a large list
    overlap_pair_sections = []
    overlap_trio_sections = []
    n_pairs = len(first_pair)
    n_trios = len(first_trio)
    for i in range(n_pairs):
        overlap_pair_sections.append([first_pair[i], second_pair[i]])
    if n_trios > 0:
        for i in range(n_trios):
            overlap_trio_sections.append([first_trio[i], second_trio[i],
                                          third_trio[i]])

    return unique_sections, overlap_pair_sections, overlap_trio_sections


# Merge overlapping sections
def merge_overlap(overlap_pair_section,
                  sdqflags=31743,
                  weight='sensitivity',
                  kind='linear'):
    '''
    Merges overlapping spectral regions. The basic workflow of this function
    is to interpolate the sections into a common wavelength table and calculate
    the weighted mean flux for each wavelength bin. If the fluxes are
    inconsistent between each other, the code can use the flux with higher SNR
    instead of the mean. If there are still outlier fluxes (compared to
    neighboring pixels), the code uses the flux from the lower SNR section
    instead.

    Parameters
    ----------
    overlap_pair_section : ``list``
        List of dictionaries containing the overlapping spectra of neighboring
        orders.

    sdqflags : ``int``, optional
        Bitwise-OR mask of serious data quality flags (``SDQFLAGS``).  Default=31743,
        which is all flags except 1024.  Flux values are averaged from non-serious-DQ
        input pixels when possible, and from serious-DQ input pixels when required.

    weight : ``str``, optional
        Defines how to merge the overlapping sections. The options currently
        implemented are ``'sensitivity'``, ``'sensitivity-dataset'`` (NET / FLUX),
        and ``'snr'`` (inverse square of the uncertainties). Default is ``'sensitivity'``.

    kind : ``str``, optional
        One of {``'linear'``, ``'zero'``, ``'slinear'``, ``'quadratic'``, ``'cubic'``}.
        Interpolation or resampling method used.
        Default is ``'linear'``.  ``'zero'``, ``'slinear'``, ``'quadratic'``
        and ``'cubic'`` refer to a spline interpolation of zeroth, first,
        second and third order, respectively.

    Returns
    -------
    overlap_merged : ``dict``
        Dictionary containing the merged overlapping spectrum.
    '''
    warnings.filterwarnings('ignore', message='Mean of empty slice', category=RuntimeWarning)  # np.nanmean

    reference = overlap_pair_section[0]
    to_shift = deepcopy(overlap_pair_section[1])

    assert np.array_equal(reference['wavelength'], sorted(reference['wavelength']))
    assert np.array_equal(to_shift['wavelength'], sorted(to_shift['wavelength']))

    # Interpolate flux and uncertainty onto the target wavelength grid:
    to_shift['flux_interpolated'] = interp1d(
        to_shift['wavelength'],
        to_shift['flux'],
        kind=kind, fill_value='extrapolate')(reference['wavelength'])
    to_shift['net_interpolated'] = interp1d(
        to_shift['wavelength'],
        to_shift['net'],
        kind=kind, fill_value='extrapolate')(reference['wavelength'])
    to_shift['uncertainty_interpolated'] = interp1d(
        to_shift['wavelength'],
        to_shift['uncertainty'],
        kind=kind, fill_value='extrapolate')(reference['wavelength'])

    to_shift['dq_regridded'] = np.zeros_like(reference['wavelength'], dtype=np.uint16)
    for i, λ_ref in enumerate(reference['wavelength']):
        if λ_ref < to_shift['wavelength'][0]:
            idx1 = idx2 = 0
        elif λ_ref > to_shift['wavelength'][-1]:
            idx1 = idx2 = len(to_shift['wavelength']) - 1
        elif λ_ref in to_shift['wavelength']:
            idx1 = idx2 = np.where(to_shift['wavelength'] == λ_ref)[0][0]
        else:
            idx1 = np.searchsorted(to_shift['wavelength'], λ_ref, side='left') - 1
            idx1 = np.clip(idx1, 0, len(to_shift['wavelength']) - 2)
            idx2 = idx1 + 1

        to_shift['dq_regridded'][i] = reduce(np.bitwise_or, to_shift['data_quality'][idx1:idx2 + 1])

    # Create NaN-masked flux arrays:
    reference['flux_with_nans'] = reference['flux'].copy()
    reference['flux_with_nans'][(reference['data_quality'] & sdqflags) != 0] = np.nan
    to_shift['flux_interpolated'][(to_shift['dq_regridded'] & sdqflags) != 0] = np.nan

    # Create NaN-masked net arrays:
    reference['net_with_nans'] = reference['net'].copy()
    reference['net_with_nans'][(reference['data_quality'] & sdqflags) != 0] = np.nan
    to_shift['net_interpolated'][(to_shift['dq_regridded'] & sdqflags) != 0] = np.nan

    # Create NaN-masked uncertainty arrays:
    reference['uncertainty_with_nans'] = reference['uncertainty'].copy()
    reference['uncertainty_with_nans']
    reference['uncertainty_with_nans'][(reference['data_quality'] & sdqflags) != 0] = np.nan
    to_shift['uncertainty_interpolated'][(to_shift['dq_regridded'] & sdqflags) != 0] = np.nan

    # Keep DQ bits associated with retained data.  If no data retained, use the bitwise-or of all:
    combined_dq = np.where(np.isnan(reference['flux_with_nans']), 0, reference['data_quality']) | \
        np.where(np.isnan(to_shift['flux_interpolated']), 0, to_shift['dq_regridded']) | \
        np.where(np.isnan(reference['flux_with_nans']) &
                 np.isnan(to_shift['flux_interpolated']),
                 reference['data_quality'] | to_shift['dq_regridded'], 0)

    # Stack arrays:
    flux_arrays = np.array([reference['flux_with_nans'], to_shift['flux_interpolated']])
    net_arrays = np.array([reference['net_with_nans'], to_shift['net_interpolated']])
    uncertainty_arrays = np.array([reference['uncertainty_with_nans'], to_shift['uncertainty_interpolated']])

    if weight == 'sensitivity':
        _, weights_ref = calculate_sensitivity(
            reference['filename'], ext=reference['ext'],
            sporder=reference['sporder'],
            output_wavelengths=reference['wavelength'])
        _, weights_to_shift = calculate_sensitivity(
            to_shift['filename'], ext=to_shift['ext'],
            sporder=to_shift['sporder'],
            output_wavelengths=reference['wavelength'])
        weights = np.array([weights_ref, weights_to_shift])
        valid_mask = ~np.isnan(flux_arrays)
        weights = np.where(valid_mask, weights, 0.)
        # Scale weights by exptime (needed if combining different observations)

        with warnings.catch_warnings(record=True) as _:  # Ignore 0 and NaN in divisors
            combined_flux = np.nansum(flux_arrays * weights, axis=0) / np.nansum(weights, axis=0)
            combined_uncertainty = np.nansum(uncertainty_arrays**2 * weights**2, axis=0)**0.5 / np.nansum(weights, axis=0)

    elif weight == 'sensitivity-dataset':
        valid_mask = ~np.isnan(flux_arrays)
        weights = np.where(valid_mask, net_arrays / flux_arrays, 0.)  # proportional to sensitivity
        # Scale weights by exptime (needed if combining different observations)

        with warnings.catch_warnings(record=True) as _:  # Ignore 0 and NaN in divisors
            combined_flux = np.nansum(flux_arrays * weights, axis=0) / np.nansum(weights, axis=0)
            combined_uncertainty = np.nansum(uncertainty_arrays**2 * weights**2, axis=0)**0.5 / np.nansum(weights, axis=0)

    elif weight in {'snr', 'inverse-variance'}:
        valid_mask = ~np.isnan(flux_arrays)
        weights = np.where(valid_mask, 1. / uncertainty_arrays**2, 0.)

        # Compute combined flux (inverse-variance weighted):
        with warnings.catch_warnings(record=True) as _:  # Ignore 0 and NaN in divisors
            combined_flux = np.nansum(flux_arrays * weights, axis=0) / np.nansum(weights, axis=0)
            combined_uncertainty = 1. / np.sqrt(np.nansum(weights, axis=0))

        # Set uncertainty to NaN where no valid data contributed:
        combined_uncertainty[np.nansum(weights, axis=0) == 0] = np.nan
        combined_flux[np.nansum(weights, axis=0) == 0] = np.nan

    elif weight == 'equal':
        combined_flux = np.nanmean([reference['flux_with_nans'], to_shift['flux_interpolated']], axis=0)
        combined_uncertainty = np.sqrt(np.nanmean(uncertainty_arrays**2, axis=0))

    elif weight == 'binary':
        combined_flux = reference['flux_with_nans']
        combined_uncertainty = reference['uncertainty']
        combined_dq = reference['data_quality']  # OVERWRITE COMBINED DQ

    else:
        raise ValueError(f'Unexpected weight value:  "{weight}"')

    overlap_merged = {
        'wavelength': reference['wavelength'],
        'flux': combined_flux,
        'uncertainty': combined_uncertainty,
        'data_quality': combined_dq, }

    return overlap_merged


# Splice the spectra
def concatenate_sections(unique_spectra_list, merged_pair_list,
                         merged_trio_list):
    """
    Concatenate the unique and the (merged) overlapping spectra.

    Parameters
    ----------
    unique_spectra_list : ``list``
        List of unique spectra.

    merged_pair_list : ``list``
        List of merged overlapping pair spectra.

    merged_trio_list : ``list``
        List of merged overlapping trio spectra.

    Returns
    -------
    spliced_wavelength : ``numpy.ndarray``
        Array containing the wavelengths in the entire spectrum.

    spliced_flux : ``numpy.ndarray``
        Array containing the fluxes in the entire spectrum.

    spliced_uncertainty : ``numpy.ndarray``
        Array containing the flux uncertainties in the entire spectrum.
    """
    n_pair_overlap = len(merged_pair_list)
    all_spectra = []

    # We always start with the first unique spectrum
    all_spectra.append(unique_spectra_list[0])
    k = 1
    for i in range(n_pair_overlap):
        all_spectra.append(merged_pair_list[i])
        try:
            all_spectra.append(merged_trio_list[i])
        except IndexError:
            all_spectra.append(unique_spectra_list[k])
            k += 1

    # There may still be some unique spectra remaining
    unique_remaining = len(unique_spectra_list) - k
    for ik in range(unique_remaining):
        all_spectra.append(unique_spectra_list[k])
        k += 1

    spliced_wavelength = \
        np.concatenate([spectrum['wavelength'] for spectrum in all_spectra])
    spliced_flux = \
        np.concatenate([spectrum['flux'] for spectrum in all_spectra])
    spliced_uncertainty = \
        np.concatenate([spectrum['uncertainty'] for spectrum in all_spectra])
    spliced_data_quality = \
        np.concatenate([spectrum['data_quality'] for spectrum in all_spectra])

    return spliced_wavelength, spliced_flux, spliced_uncertainty, \
        spliced_data_quality


# The splice pipeline does everything
def splice(x1d_input, ext=1, update_fits=False, output_file=None, weight='sensitivity',
           kind='linear', sdqflags=31743,
           truncate_edge_left=5, truncate_edge_right=5):
    """
    The main workhorse of the package. This pipeline performs all the steps
    necessary to merge overlapping spectral sections and splice them with the
    unique sections.

    Parameters
    ----------
    x1d_input : ``str``
        Path and name of the ``*_x1d.fits`` file containing the spectrum.

    ext : ``int``, optional
        FITS extension to read from the X1D file.  Most X1D files contain only
        one SCI extension (i.e., ``REPEATOBS=1``).  Default is ``1``.

    update_fits : ``bool``, optional
        Use carefully, since it can modify fits files permanently. Parameter
        that decides whether to update the ``*_x1d.fits`` file with a new
        extension containing the spliced spectrum.

    output_file : ``str`` or ``None``, optional
        String containing the location to save the output spectrum as an ascii
        file. If ``None``, no output file is saved and the code returns an
        Astropy Table instead. Default is ``None``.

    weight : ``str``, optional
        Defines how to merge the overlapping sections. The options currently
        implemented are ``'sensitivity'`` and ``'snr'`` (inverse square of the
        uncertainties). Default is ``'sensitivity'``.

    kind : ``str``, optional
        One of {``'linear'``, ``'zero'``, ``'slinear'``, ``'quadratic'``, ``'cubic'``}.
        Interpolation or resampling method used.
        Default is ``'linear'``.  ``'zero'``, ``'slinear'``, ``'quadratic'``
        and ``'cubic'`` refer to a spline interpolation of zeroth, first,
        second and third order, respectively.

    sdqflags : ``int``, optional
        Serious data quality flags.  This is a bitwise-or of the flag values
        to use when excluding data from its source.

    truncate_edge_left : ``int``, optional
        Set the number of low-resolution pixels at the left edge of the detector
        where the spectra should be truncated. If ``None``, then no truncation
        is applied. Default is ``5``.

    truncate_edge_right : ``int``, optional
        Set the number of low-resolution pixels at the right edge of the
        detector where the spectra should be truncated. If ``0``, then no
        truncation is applied. Default is ``5``.

    Returns
    -------
    spliced_spectrum_table : ``astropy.Table`` object
        Astropy Table containing the spliced spectrum. Only returned if
        ``output_file`` is ``None``.
    """
    # Read the data
    sections = read_spectrum(x1d_input, ext=ext, truncate_edge_left=truncate_edge_left,
                             truncate_edge_right=truncate_edge_right)

    unique_sections, overlap_pair_sections, overlap_trio_sections = \
        find_overlap(sections, extent=1024 - truncate_edge_left - truncate_edge_right)

    # Merge the overlapping spectral sections
    merged_pairs = [
        merge_overlap(overlap_pair_sections[k], sdqflags, weight, kind=kind)
        for k in range(len(overlap_pair_sections))
    ]

    if len(overlap_trio_sections) > 0:
        merged_trios = [
            merge_overlap(overlap_trio_sections[k], sdqflags, weight, kind=kind)
            for k in range(len(overlap_trio_sections))
        ]
    else:
        merged_trios = []

    # By now we have three lists: unique_sections, merged_pairs and
    # merged_trios. The next step is to concatenate everything in the correct
    # order.

    # Finally splice the unique and merged sections
    wavelength, flux, uncertainty, dq = concatenate_sections(unique_sections,
                                                             merged_pairs,
                                                             merged_trios)

    # Instantiate the spectrum dictionary
    spectrum_dict = \
        {'WAVELENGTH': wavelength, 'FLUX': flux, 'ERROR': uncertainty, 'DQ': dq}
    spliced_spectrum_table = Table(spectrum_dict)
    table_hdu = fits.BinTableHDU(spliced_spectrum_table)

    # This feature modifies the fits input file! Use carefully!
    if update_fits:
        with fits.open(x1d_input, mode='update') as hdul:
            hdul.append(table_hdu)

    # Return or output the result
    if output_file is None:
        return spliced_spectrum_table
    else:
        spliced_spectrum_table.write(output_file, format='fits', overwrite=True)
