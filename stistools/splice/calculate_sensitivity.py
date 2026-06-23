import warnings
import numpy as np
from astropy.io import fits
from astropy.table import Table
from stistools.r_util import expandFileName


__doc__ = r"""
The ``calculate_sensitivity`` module is intended to estimate the instrument sensitivity
to facilitate weighting in ``splice``.  While sensitivity is proportional to NET / FLUX,
direct evaluation of this ratio does not work for pixels containing zero counts.  Thus,
we directly evaluate the sensitivity using information from STIS reference files in the
header of X1D/SX1 files, including the time- and temperature-dependence, blaze shift,
aperture throughput, extraction height correction, etc.

While throughputs are available via the synphot/stsynphot packages, these data do not
provide accurate inter-order echelle throughput levels.  Furthermore, users would need
to download and install synphot data.

Notes
-----
The results of ``calculate_sensitivity`` are not used by ``calstis`` to directly
flux-calibrate data.  We do not recommend its use for this purpose.  The algorithms
here attempt to recreate the logic found in the ``hstcal`` package.
"""


SEC_PER_DAY = 24 * 60 * 60
HC = 1.9864e-8  # Plank's constant * speed of light (erg * Å)
A_HST = 45_238.93416  # effective area of HST primary mirror (cm^2)


def eval_tds(dataset, ext=1, tempcorr=True, output_wavelengths=None):
    '''Evaluate the Time-Dependent Sensitivity (TDS) at a specified observation's MJD
    and temperature, and then optionally interpolate to a specified wavelength grid (in Å).

    Parameters
    ----------
    dataset : str
        Path to a STIS x1d or sx1 FITS file

    ext : int, optional
        dataset extension from which to get the temperature.

    tempcorr : bool
        If True (default), apply the temperature-sensitivity correction

    output_wavelengths : np.ndarray, optional
        Wavelength array (Å) onto which the TDS correction
        should be interpolated. If omitted, the native TDSTAB
        wavelength grid is returned.

    Returns
    -------
    wavelengths : np.ndarray
        Wavelength array in Å.

    throughputs : np.ndarray
        TDS throughput correction evaluated at `wavelengths`.
    '''
    with fits.open(dataset) as f:
        mjd = f[0].header['TEXPSTRT']

        detector = f[0].header['DETECTOR']
        opt_elem = f[0].header['OPT_ELEM']
        if tempcorr:
            if detector == 'FUV-MAMA':
                temperature = f[ext].header['OM1CAT']
            elif detector == 'NUV-MAMA':
                temperature = f[ext].header['OM2CAT']
            elif detector == 'CCD':
                raise NotImplementedError('Does not yet support CCD')
                # check if on side-2 electronics...
            else:
                raise ValueError(f'Unexpected detector="{detector}"')
        else:
            temperature = np.nan

        tdstab_filename = expandFileName(f[0].header['TDSTAB'])

    tds = Table.read(tdstab_filename, hdu=1, unit_parse_strict='silent')
    tds['OPT_ELEM'] = [x.strip() for x in tds['OPT_ELEM']]

    try:
        tds = tds[tds['OPT_ELEM'] == opt_elem][0]
    except IndexError as e:
        raise ValueError(f'Unexpected OPT_ELEM="{opt_elem}" for TDSTAB="{tdstab_filename}"') from e

    slope = np.asarray(tds['SLOPE'][:tds['NT'], :tds['NWL']], dtype=float)
    tds_wavelength = tds['WAVELENGTH'][:tds['NWL']]

    dates = tds['TIME'][:tds['NT']]
    # Add a future date beyond the observation date to enable interpolation:
    with warnings.catch_warnings(record=True) as _:
        future = max([tds['TIME'][:tds['NT']][-1] + 365,
                      mjd + 365])
    dates = np.asarray(np.append(dates, [future]), dtype=float)

    # Replicate last element:
    slope = np.resize(slope, (tds['NT'] + 1, tds['NWL']))
    slope[-1, :] = slope[-2, :]

    y = np.zeros((tds['NT'] + 1, tds['NWL']), dtype=float)
    for t_index, t in enumerate(dates):
        factor = np.ones((tds['NWL']), dtype=float)
        for i in range(tds['NWL']):
            for j in range(tds['NT']):
                slope_j = slope[j, i] / 365.25 / 100.0  # %/yr
                if (j == tds['NT'] - 1) or (t <= tds['TIME'][j + 1]):
                    factor[i] += (t - tds['TIME'][j]) * slope_j
                    break
                else:
                    factor[i] += (tds['TIME'][j + 1] - tds['TIME'][j]) * slope_j
        y[t_index, :] = factor

    if tempcorr and not np.isnan(temperature):
        y *= (1. + tds['TEMPSENS'][:tds['NWL']] * (temperature - tds['REFTEMP']))

    # Interpolate each wavelength bin to the observation date:
    yy = np.zeros((len(tds_wavelength)), dtype=float)
    for i, wl in enumerate(tds_wavelength):
        yy[i] = np.interp(mjd, dates, y[:, i])

    # Interpolate the TDS factor to the desired wavelength grid:
    if output_wavelengths is not None:
        yy = np.interp(output_wavelengths, tds_wavelength, yy)
    else:
        output_wavelengths = tds_wavelength

    return output_wavelengths, yy


def aperture_throughput(dataset, output_wavelengths=None):
    '''Locate, parse, and interpolate the aperture throughput correction for a
    STIS dataset.

    Parameters
    ----------
    dataset : str
        Path to a STIS x1d or sx1 FITS file

    output_wavelengths : np.ndarray, optional
        Wavelength array (Å) onto which the throughput correction
        should be interpolated. If omitted, the native APERTAB
        wavelength grid is returned.

    Returns
    -------
    wavelengths : np.ndarray
        Wavelength array in Å.

    throughputs : np.ndarray
        Aperture throughput correction evaluated at `wavelengths`.
    '''
    dataset_aperture = fits.getval(dataset, ext=0, keyword='APERTURE')
    apertab_filename = expandFileName(fits.getval(dataset, ext=0, keyword='APERTAB'))
    apertab = Table.read(apertab_filename, hdu=1, unit_parse_strict='silent')

    # Filter APERTAB to matching rows:
    apertab = apertab[[x['APERTURE'].strip() == dataset_aperture for x in apertab]]
    if len(apertab) != 1:
        raise ValueError(f'Unique APERTAB entry not found for aperture="{dataset_aperture}"')
    apertab = apertab[0]

    wavelength = apertab['WAVELENGTH'][:apertab['NELEM']]
    throughput = apertab['THROUGHPUT'][:apertab['NELEM']]

    if output_wavelengths is not None:
        throughput = np.interp(output_wavelengths, wavelength, throughput)
        wavelength = output_wavelengths

    return wavelength, throughput


def extraction_height_correction(dataset, ext=1, output_wavelengths=None):
    '''Locate, parse, and interpolate the PCTAB extraction height correction for a
    STIS dataset.

    Parameters
    ----------
    dataset : str
        Path to a STIS x1d or sx1 FITS file

    ext : int
        Extension of dataset.  Used to look up extraction height used.  Default=1.

    output_wavelengths : np.ndarray, optional
        Wavelength array (Å) onto which the throughput correction
        should be interpolated. If omitted, the native wavelength
        grid is returned.

    Returns
    -------
    wavelengths : np.ndarray
        Wavelength array in Å.

    throughputs : np.ndarray
        Extraction throughput correction evaluated at `wavelengths`.
    '''
    with fits.open(dataset) as f:
        pctab_filename = expandFileName(fits.getval(dataset, ext=0, keyword='PCTAB'))
        pctab = Table.read(pctab_filename, hdu=1, unit_parse_strict='silent')
        pctab = pctab[[x['CENWAVE'] in {-1, f[0].header['CENWAVE']} and
                      (x['APERTURE'] == f[0].header['APERTURE']) for x in pctab]]
        pctab.sort('EXTRHEIGHT')

        pctab_inf = pctab[pctab['EXTRHEIGHT'] == 600][0]
        wavelength_inf = pctab_inf['WAVELENGTH'][:pctab_inf['NELEM']]
        throughput_inf = pctab_inf['THROUGHPUT'][:pctab_inf['NELEM']]

        try:
            extrsize = set(f[ext].data['EXTRSIZE'])
        except IndexError as e:
            raise ValueError(f"Dataset {dataset} has no ext={ext}") from e
        if len(extrsize) > 1:
            raise ValueError('Multiple EXTRSIZEs present')
        extrsize = float(list(extrsize)[0])
        if extrsize in pctab['EXTRHEIGHT']:
            pctab_match = pctab[pctab['EXTRHEIGHT'] == extrsize][0]
            wavelength = pctab_match['WAVELENGTH'][:pctab_match['NELEM']]
            throughput = pctab_match['THROUGHPUT'][:pctab_match['NELEM']]
        elif (extrsize >= pctab['EXTRHEIGHT'].min()) and (extrsize <= pctab['EXTRHEIGHT'].max()):
            extrsize_below = pctab['EXTRHEIGHT'][pctab['EXTRHEIGHT'] < extrsize].max()
            extrsize_above = pctab['EXTRHEIGHT'][pctab['EXTRHEIGHT'] >= extrsize].min()
            pctab_match_below = pctab[pctab['EXTRHEIGHT'] == extrsize_below][0]
            pctab_match_above = pctab[pctab['EXTRHEIGHT'] == extrsize_above][0]
            wavelength_below = pctab_match_below['WAVELENGTH'][:pctab_match_below['NELEM']]
            wavelength_above = pctab_match_above['WAVELENGTH'][:pctab_match_above['NELEM']]
            throughput_below = pctab_match_below['THROUGHPUT'][:pctab_match_below['NELEM']]
            throughput_above = pctab_match_above['THROUGHPUT'][:pctab_match_above['NELEM']]
            if (len(wavelength_below) != len(wavelength_above)) and (wavelength_below != wavelength_above):
                raise ValueError(f'PCTAB has different wavelengths for EXTRHEIGHTs {extrsize_below} and {extrsize_above}.')
            wavelength = wavelength_below
            # Interpolate throughput at each wavelength bin:
            throughput = np.zeros((len(wavelength_below)), dtype=float)
            for i, (thr_below, thr_above) in enumerate(zip(wavelength_below, throughput_below, throughput_above)):
                throughput[i] = np.interp(extrsize, [extrsize_below, extrsize_above], [thr_below, thr_above])
        else:
            raise ValueError(f'EXTRSIZE={extrsize} is beyond the extent of the PCTAB')

    if (len(wavelength) != len(wavelength_inf)) and (wavelength != wavelength_inf):
        raise ValueError(f'PCTAB has different wavelengths for EXTRHEIGHTs {extrsize} and 600.')
    throughput_ratio = throughput_inf / throughput

    if output_wavelengths is not None:
        throughput_ratio = np.interp(output_wavelengths, wavelength, throughput_ratio)
        wavelength = output_wavelengths

    return wavelength, throughput_ratio


def calculate_sensitivity(dataset, ext, sporder, blazecorr=True, output_wavelengths=None):
    '''Calculate the instrument sensitivity corresponding to a specified STIS X1D/SX1
    dataset.

    Parameters
    ----------
    dataset : str
        Path to a STIS x1d or sx1 FITS file

    ext : int
        Extension of dataset

    sporder : int
        Spectral order number to be calculated

    blazecorr : bool
        If True (default), perform the blaze shift correction to the throughput.

    output_wavelengths : np.ndarray, optional
        Wavelength array (Å) onto which the throughput should be
        interpolated. If omitted, the PHOTTAB's native wavelength
        grid is returned.

    Returns
    -------
    wavelengths : np.ndarray
        Wavelength array in Å.

    sensitivity : np.ndarray
        Sensitivity evaluated at `wavelengths`.

    Notes
    -----
    Resulting sensitivity matches NET/FLUX well, but is not subject to undefined behavior
    when the count rate is 0 for a pixel.
    '''
    with fits.open(dataset) as f:
        gain = f[0].header.get('ATODGAIN', 1.)  # 1 for MAMAs
        aperture = f[0].header['APERTURE']
        opt_elem = f[0].header['OPT_ELEM']
        cenwave = f[0].header['CENWAVE']

        phottab_filename = expandFileName(f[0].header['PHOTTAB'])
        phottab = Table.read(phottab_filename, hdu=1, unit_parse_strict='silent')
        # Filter PHOTTAB to matching rows:
        phottab = phottab[(phottab['OPT_ELEM'] == opt_elem) & (phottab['CENWAVE'] == cenwave)]

        if f[ext].header.get('DOPPON', False):
            doppmag = f[ext].header['DOPPMAG']
            doppzero = f[ext].header['DOPPZERO']
            orbitper = f[ext].header['ORBITPER'] / SEC_PER_DAY
            expstart = f[ext].header['EXPSTART']
            expend = f[ext].header['EXPEND']
            if (expend - expstart) <= 1e-4:  # 8.64 s
                mid_time = (expend + expstart) / 2.
                doppler_shift = doppmag * np.sin((mid_time - doppzero) * 2. * np.pi / orbitper)
            else:
                doppler_shift = doppmag * orbitper / (2. * np.pi) * \
                    (np.cos((expstart - doppzero) * 2. * np.pi / orbitper) -
                     np.cos((expend - doppzero) * 2. * np.pi / orbitper)) \
                    / (expend - expstart)

            doppler_shift /= 2  # convert to low-res pixels
        else:
            doppler_shift = 0.

        v_helio = f[ext].header.get('V_HELIO', 0.)  # km/s
        heliocentric_factor = 1. + v_helio / 299_792.458  # 1 + v/c

        order = phottab[phottab['SPORDER'] == sporder][0]

        x1d_order = f[ext].data[f[ext].data['SPORDER'] == order['SPORDER']][0]
        # Convert the x1d file's wavelengths back to the detector frame:
        x1d_λ = x1d_order['WAVELENGTH'] * heliocentric_factor
        disp = float(np.diff(x1d_λ[511:514]).mean())
        x1d_λ += doppler_shift * disp

        if blazecorr:
            # Determine values for the reference order:
            ref_order = f[ext].data[f[ext].data['SPORDER'] == order['REFORD']][0]
            # Convert the reference order's wavelengths back to the detector frame:
            ref_λ = ref_order['WAVELENGTH'] * heliocentric_factor  # Remove HELCORR
            disp_ref = float(np.diff(ref_λ[511:514]).mean())
            ref_λ += doppler_shift * disp_ref  # Remove DOPPCORR
            # reference order's wavelength and Y-position (1-indexed) at x-center:
            obsw = float(ref_λ[512])
            obsy = float(ref_order['A2CENTER'])

            # Recalculate disp to account for the small heliocentric_factor:
            disp = float(np.diff(ref_order['WAVELENGTH'][511:514]).mean())
            Δx = (order['REFWAV'] - obsw) / disp
            Δy = obsy - order['REFY']
            Δt = expstart - order['REFMJD']
            # Calculate the blaze shift offset:
            bshift = order['BSHIFT_VS_X'] * Δx + order['BSHIFT_VS_Y'] * Δy + order['BSHIFT_VS_T'] * Δt + order['BSHIFT_OFFSET']

            # Offset wavelengths by the blaze shift:
            wavephtshift_fit = order['WAVELENGTH'] + bshift * disp * ref_order['SPORDER'] / order['SPORDER']
        else:
            wavephtshift_fit = order['WAVELENGTH']

        # Interpolate the PHOTTAB's throughput onto the x1d file's wavelength grid,
        # shifting by the blaze shift:
        throughput = np.interp(x1d_λ, wavephtshift_fit, order['THROUGHPUT'])

        # Calculate correction factors
        _, tds = eval_tds(dataset, ext=ext, output_wavelengths=x1d_λ)  # time- and temperature-dependent sensitivities
        _, H = extraction_height_correction(dataset, ext=ext, output_wavelengths=x1d_λ)
        wavelength, aper = aperture_throughput(dataset, output_wavelengths=x1d_λ)

        Δλ = np.diff(x1d_λ)
        Δλ = np.concatenate([Δλ, [Δλ[-1]]])

        sensitivity = throughput * A_HST * x1d_λ * Δλ * tds * aper / (HC * gain * H)

        if output_wavelengths is not None:
            sensitivity = np.interp(output_wavelengths, wavelength, sensitivity)
            wavelength = output_wavelengths

    return wavelength, sensitivity
