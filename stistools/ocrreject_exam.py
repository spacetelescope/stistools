#!/usr/bin/env python

import os
import warnings
import argparse

import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm as colormap
from matplotlib import colors
from matplotlib import pyplot as plt
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

__doc__ = """
    Checks STIS CCD 1D spectroscopic data for cosmic ray overflagging.

    Examples
    --------

    In Python:

    >>> import stistools
    >>> stistools.ocrreject_exam.ocrreject_exam("odvkl1040", plot=True)

    .. code-block:: python

       [{'rootname': 'odvkl1040',
       'extr_fracs': array([0.31530762, 0.32006836]),
       'outside_fracs': array([0.00884673, 0.00810278]),
       'ratios': array([35.64113429, 39.50106762]),
       'avg_extr_frac': 0.31768798828125,
       'avg_outside_frac': 0.008474755474901575,
       'avg_ratio': 37.486389928547126}]

    .. image:: odvkl1040_stacked.png
      :width: 600
      :alt: Stacked example ocrreject_exam plot output

    |

    .. image:: odvkl1040_splits.png
      :width: 600
      :alt: Split example ocrreject_exam plot output

    From command line:

    .. code-block:: none

       ocrreject_exam -h
       usage: ocrreject_exam [-h] [-d DATA_DIR] [-p] [-o PLOT_DIR] [-i] obs_id [obs_id ...]

       Calculate fractions of cosmic ray rejected pixels inside and outside of an extraction box to test for CR algorithm failures.

       positional arguments:
       obs_id       observation id(s) in ipppssoot format

       options:
       -h, --help   show this help message and exit
       -d DATA_DIR  directory containing observation flt and sx1/x1d files. Defaults to current working directory.
       -p           option to create diagnostic plots
       -o PLOT_DIR  output directory to store diagnostic plots if plot=True. Defaults to data_dir.
       -i           option to create zoomable html plots instead of static pngs. Defaults to False and requires Plotly if True

       v1.0; Written by Matt Dallas, Joleen Carlberg, Sean Lockwood, STScI, December 2024.
    """

__taskname__ = "ocrreject_exam"
__version__  = "1.0"
__vdate__    = "09-December-2024"
__author__   = "Matt Dallas, Joleen Carlberg, Sean Lockwood, STScI, December 2024."


class BoxExtended(Exception):
    def __init__(self, message='Extraction box extends beyond frame'):
        super().__init__(message)


def ocrreject_exam(obs_ids, data_dir='.', plot=False, plot_dir=None, interactive=False,
                   verbose=False):
    """Compares the rate of cosmic rays in the extraction box and everywhere else 
    in a CCD spectroscopic image. Based on crrej_exam from `STIS ISR 2019-02 
    <https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/stis/documentation/instrument-science-reports/_documents/201902.pdf>`_.

    Higher ratios of cosmic ray rates in the extraction box to the rest of the image
    may indicate the need to rerun stistools.ocrreject() with different parameters.

    Parameters
    ----------
    obs_ids: iter of str or str
        One of more STIS observation ID rootnames in ipppssoot format (e.g. odvkl1040). 

    data_dir: str
        Directory containing both the flat fielded (_flt.fits) and extracted 
        spectrum (_sx1.fits or _x1d.fits) files of the observation if using obs_ids argument. 
        Defaults to current working directory.

    plot: bool
        Option to generate diagnostic plots, default=False

    plot_dir: str or None
        Directory to save diagnostic plots in if plot=True. 
        Defaults to data_dir parameter

    interactive: bool
        Option to generate zoomable html plots using plotly, default=False

    verbose: bool
        Option to print some results

    Returns
    -------
    results: list of dict

        - ``rootname``: obs_id
        - ``extr_fracs``: cosmic ray rejection rates in the extraction boxes for each CR-SPLIT
        - ``outside_fracs``: cosmic ray rejection rates outside the extraction boxes for each CR-SPLIT
        - ``ratios``: ``extr_fracs``/``outside_fracs``
        - ``avg_extr_frac``: The average of ``extr_fracs``
        - ``avg_outside_frac``: The average of ``outside_fracs``
        - ``avg_ratio``: ``avg_extr_frac``/``avg_outside_frac``
        
    If called from the command line, prints the avg extraction, outside, and ratio values for quick verification.
    """
    if isinstance(obs_ids, (str,)):
        obs_ids = [obs_ids]

    result_list = []
    for obs_id in obs_ids:
        flt_file = os.path.join(data_dir, obs_id.lower() + '_flt.fits')
        if not os.access(flt_file, os.R_OK):
            raise FileNotFoundError(f"FLT file for {obs_id} not found in '{data_dir}'")

        sx1_file = os.path.join(data_dir, obs_id.lower() + '_sx1.fits')
        x1d_file = os.path.join(data_dir, obs_id.lower() + '_x1d.fits')
        if os.access(sx1_file, os.F_OK) and os.access(x1d_file, os.F_OK):
            warnings.warn(f"Found both SX1 and X1D files for {obs_id} in '{data_dir}', defaulting to use SX1")

        if not os.access(sx1_file, os.R_OK):
            sx1_file = x1d_file
            if not os.access(sx1_file, os.R_OK):
                raise FileNotFoundError(f"SX1/X1D file for {obs_id} not found in '{data_dir}'")

        if plot and plot_dir is None:
            plot_dir = data_dir

        with fits.open(flt_file) as flt_hdul:
            # When first opening the flt file check that it is even of a CCD exposure
            try:
                instrument = flt_hdul[0].header['INSTRUME']
                detector   = flt_hdul[0].header['DETECTOR']
                obsmode    = flt_hdul[0].header['OBSMODE']
                nextend    = flt_hdul[0].header['NEXTEND']

                if (instrument.strip() != 'STIS') or \
                   (detector.strip() != 'CCD')    or \
                   (obsmode.strip() != 'ACCUM')   or \
                   (nextend < 6):
                    raise ValueError
            except (KeyError, ValueError,) as e:
                raise ValueError(f"{flt_file}: Not a STIS/CCD/ACCUM file with ≥2 SCI extensions.") from e

            # If it is, proceed to check that the number of sci extensions matches the number of crsplits
            propid = flt_hdul[0].header['PROPOSID']
            rootname = flt_hdul[0].header['ROOTNAME']
            nrptexp_num = flt_hdul[0].header['NRPTEXP']
            crsplit_num = flt_hdul[0].header['CRSPLIT']
            sci_num = len([hdu.name for hdu in flt_hdul if "SCI" in hdu.name]) # Counts the number of sci extensions

        if (crsplit_num * nrptexp_num - sci_num) != 0:
            raise ValueError(f"cr-split or nrptexp value in flt header does not match the number of sci extentsions for {obs_id}")

        # If all checks above passed, calculate cr fraction in and out of the extraction box
        spec = fits.getdata(sx1_file, ext=1)[0]

        extrlocy = spec['EXTRLOCY'] - 1 # y coords of the middle of the extraction box
        del_pix = spec['EXTRSIZE'] / 2. # value the extraction box extends above or below extrlocy
        box_upper = np.ceil(extrlocy + del_pix).astype(int) # Ints of pixel values above end of trace bc python is upper bound exclusive
        box_lower = np.floor(extrlocy - del_pix).astype(int) # Ints of pixel values below end of trace

        # Fill each of these lists with values for each cr split
        extr_fracs = [] # fraction of pixels flagged as cr inside the extraction box for each split
        outside_fracs = [] # fraction of pixels flagged as cr outside the extraction box for each split
        cr_rejected_locs = [] # 2d array of 1s where a cr exists and 0 elsewhere
        exposure_times = []

        with fits.open(flt_file) as flt_hdul:
            flt_shape = flt_hdul['sci', 1].data.shape # shape of the data

            # Check that the extraction box doesn't extend beyond the image: this breaks the method
            if np.any(box_lower < 0) or np.any(box_upper - 1 > flt_shape[0]): # Subtract 1 because the box extends to the value of the pixel before
                raise BoxExtended(f"Extraction box coords extend above or below the cosmic ray subexposures for {propid}")

            extr_mask = np.zeros(flt_shape)
            outside_mask = np.ones(flt_shape)

            for column in range(0, flt_shape[1]):
                extr_mask[box_lower[column]:box_upper[column], column] = 1 # 1s inside the extraction box, 0s outside
                outside_mask[box_lower[column]:box_upper[column], column] = 0 # 0s inside the extraction box, 1s outside

            n_extr = np.count_nonzero(extr_mask) # number of pixels inside the extraction box
            n_outside = np.count_nonzero(outside_mask) # number of pixels outside the extraction box

            for i, hdu in enumerate(flt_hdul):
                if hdu.name == 'SCI':
                    exposure_times.append(hdu.header['EXPTIME'])
                    dq_array = flt_hdul[i + 2].data # dq array corresponding to each sci extentsion

                    extr_rej_pix = np.zeros(flt_shape) # 2d array where there is a 1 if a pixel inside the extraction box is marked as a cr
                    np.place(extr_rej_pix, (extr_mask == 1) & (dq_array & 2**13 != 0), 1)

                    outside_rej_pix = np.zeros(flt_shape) # 2d array where there is a 1 if a pixel outside the extraction box is marked as a cr
                    np.place(outside_rej_pix, (outside_mask == 1) & (dq_array & 2**13 != 0), 1)

                    extr_cr_count = np.count_nonzero(extr_rej_pix)
                    outside_cr_count = np.count_nonzero(outside_rej_pix)

                    extr_fracs.append(extr_cr_count / n_extr)
                    outside_fracs.append(outside_cr_count / n_outside)

                    cr_rejected_pix = extr_rej_pix + outside_rej_pix
                    cr_rejected_locs.append(cr_rejected_pix)

        extr_fracs = np.asarray(extr_fracs)
        outside_fracs = np.asarray(outside_fracs)
        ratios = extr_fracs / outside_fracs # ratio of extraction to outside the box in each image

        avg_extr_frac = float(np.sum(extr_fracs) / len(extr_fracs)) # Average fraction of crs inside extraction box
        avg_outside_frac = float(np.sum(outside_fracs) / len(outside_fracs)) # Average fraction of crs outside extraction box
        avg_ratio = float(avg_extr_frac / avg_outside_frac) # Average ratio of the stack

        results = {
            'rootname'         : obs_id,
            'extr_fracs'       : extr_fracs,
            'outside_fracs'    : outside_fracs,
            'ratios'           : ratios,
            'avg_extr_frac'    : avg_extr_frac,
            'avg_outside_frac' : avg_outside_frac,
            'avg_ratio'        : avg_ratio,}

        if plot and (not interactive or not HAS_PLOTLY): # case with interactive == False
            if not HAS_PLOTLY and interactive:
                warnings.warn('Plotly required for interactive plotting, using matplotlib and static pngs.')
                interactive = False

            cr_rejected_stack = np.sum(cr_rejected_locs, axis=0) # stack all located crs on top of eachother
            stacked_exposure_time = sum(exposure_times)
            stack_plot(cr_rejected_stack, box_lower, box_upper, len(cr_rejected_locs), stacked_exposure_time,
                rootname, propid, plot_dir, interactive=interactive)
            split_plot(cr_rejected_locs, box_lower, box_upper, len(cr_rejected_locs), exposure_times,
                stacked_exposure_time, rootname, propid, plot_dir, interactive=interactive)

        elif plot and interactive and HAS_PLOTLY: # case with interactive == True and Plotly is installed
            cr_rejected_stack = np.sum(cr_rejected_locs, axis=0) # stack all located crs on top of each other
            stacked_exposure_time = sum(exposure_times)
            stack_plot(cr_rejected_stack, box_lower, box_upper, len(cr_rejected_locs), stacked_exposure_time,
                       rootname, propid, plot_dir, interactive=interactive)
            split_plot(cr_rejected_locs, box_lower, box_upper, len(cr_rejected_locs), exposure_times,
                stacked_exposure_time, rootname, propid, plot_dir, interactive=interactive)

        if verbose:
            print(f"\nFor {obs_id}")
            print(f"Average across all extraction boxes: {results['avg_extr_frac']:.1%}")
            print(f"Average across all external regions: {results['avg_outside_frac']:.1%}")
            print(f"Average ratio between the two: {results['avg_ratio']:.2f}")

        result_list.append(results)

    return result_list


# Plotting-specific functions:
def _gen_color(cmap, n):
    """Generates n distinct colors from a given colormap. 

    Based on mycolorpy's gen_color() from https://github.com/binodbhttr/mycolorpy
    """
    colorlist = []

    for c in cmap.colors[0:n]:
        clr = colors.rgb2hex(c) # convert to hex
        colorlist.append(str(clr)) # create a list of these colors

    colorlist.pop(0) # Make it light grey rather than black at the beginning (I think it's easier on the eyes)
    colorlist.insert(0, '#F5F5F5')

    return colorlist


def _discrete_colorscale(bvals, colors):
    """Takes desired boundary values and colors from a matplotlib colorplot and makes a Plotly colorscale.

    Based on discrete_colorscale() from https://community.plotly.com/t/colors-for-discrete-ranges-in-heatmaps/7780
    """
    if len(bvals) != len(colors) + 1:
        raise ValueError('len(boundary values) should be equal to len(colors) + 1')
    bvals = sorted(bvals)
    nvals = [(v - bvals[0]) / (bvals[-1] - bvals[0]) for v in bvals]  # normalized values

    dcolorscale = [] # discrete colorscale
    for k, color in enumerate(colors):
        dcolorscale.extend([[nvals[k], color], [nvals[k+1], color]])

    return dcolorscale


def _generate_intervals(n, divisions):
    """Creates a list of strings that are the positions requred for centering an evenly spaced colorbar in Plotly
    """
    result = np.linspace(0, n, divisions, endpoint=False)
    offset = (result[1] - result[0]) / 2
    result = result + offset
    result = [str(x) for x in list(result)[:len(list(result))]]

    return result


def stack_plot(stack_image, box_lower, box_upper, split_num, texpt, obs_id, propid, plot_dir,
               interactive):
    """Creates a visualization of where CR pixels are in a stacked image

    Parameters
    ----------
    stack_image: array
        2d array to plot

    box_lower: array
        1d array of ints of the bottom of the extraction box (0 indexed)

    box_upper: array
        1d array of ints of the top of the extraction box (0 indexed)

    split_num: int
        Number of splits in the stack

    texpt: float
        Value of total exposure time

    obs_id: str
        ipppssoot of observation

    propid: int
        proposal id of observation

    plot_dir: str
        Directory to save plot in

    interactive: bool 
        If True, uses Plotly to create an interactive zoomable html plot
    """
    stack_shape = stack_image.shape
    max_stack_value = int(np.max(stack_image)) # This is usually equal to stack_shape,
    # in the case where a cr pixel is not in all splits at the same location this value should be used

    color_list = [
        'k', 'tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:cyan', 'tab:olive', 'tab:purple',
        'tab:pink', 'tab:brown', 'tab:grey', 'darkkhaki', 'gold', 'lightskyblue', 'peru', 'slateblue',
        'darkolivegreen', 'mediumseagreen', 'tomato', 'paleturquoise', 'lightgreen', 'chocolate',
        'yellowgreen', 'darksalmon', 'olive', 'darkgoldenrod', 'firebrick', 'teal', 'magenta',
        'mediumaquamarine', 'darkslategrey', 'blueviolet', 'peachpuff',]
    # hard-coded to 32 values, this should cover all cr split numbers

    custom_cmap = colors.ListedColormap(color_list)
    cmap = colors.ListedColormap(_gen_color(custom_cmap, max_stack_value + 1))
    bounds = np.arange(max_stack_value + 2)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    if not interactive:
        # create matplotlib image
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(9, 20*(9/41)),
            gridspec_kw={'width_ratios': [1, 1, 0.05], 'height_ratios': [1]})

        for axis in [ax1, ax2]:
            axis.imshow(stack_image, interpolation='none', origin='lower',
                extent=(0, stack_shape[1], 0, stack_shape[0]), cmap=cmap, norm=norm, aspect='auto')
            axis.step(np.arange(len(box_upper)), box_upper, color='#222222', where='post', lw=0.7, alpha=0.7, ls='--')
            axis.step(np.arange(len(box_lower)), box_lower, color='#222222', where='post', lw=0.7, alpha=0.7, ls='--')

        ax1.set_title('Full image')

        # If it is a large enough image, zoom the 2nd subplot around the extraction box region
        if ((stack_shape[0] - max(box_upper)) > 20) and (min(box_lower) > 20):
            ax2.set_ylim([(min(box_lower)-20),(max(box_upper)+20)])
            ax2.set_title('zoomed to 20 pixels above/below extraction box')

        # Otherwise just don't zoom in at all
        else:
            ax2.set_title('full image already 20 pixels above/below extraction box')

        cb = fig.colorbar(colormap.ScalarMappable(norm=norm, cmap=cmap), cax=ax3,
            label='# times flagged as CR', ticks=np.arange(max_stack_value, max_stack_value + 2) - 0.5)
        cb.set_ticklabels(np.arange(max_stack_value, max_stack_value+2)-1)

        fig.suptitle(f"CR flagged pixels in stacked image: {obs_id}\n Proposal {propid!s}, " \
                     f"exposure time {texpt:.2f}, {split_num!s} subexposures")
        fig.tight_layout()

        plot_name = obs_id + '_stacked.png'
        file_path = os.path.join(plot_dir, plot_name)
        plt.savefig(file_path, dpi=150, bbox_inches='tight')
        plt.close()

    else:
        # Create Plotly image
        fig = go.Figure()

        # calculate required x and y range, colorbar info, and figure titles
        x = np.arange(start=0, stop=stack_shape[1] + 1, step=1)
        y = np.arange(start=0, stop=stack_shape[0] + 1, step=1)

        dcolorsc = _discrete_colorscale(bvals=list(bounds), colors=cmap.colors)

        ticktext = [str(x) for x in list(bounds)[:len(list(bounds)) - 1]]
        tickvals = _generate_intervals(len(ticktext) - 1, len(ticktext))

        title_text = f"CR flagged pixels in stacked image: {obs_id}<br>Proposal {propid!s}, " \
                     f"exposure time {texpt:.2f}, {split_num!s} subexposures"
        plot_name = obs_id + '_stacked.html'
        file_path = os.path.join(plot_dir, plot_name)

        # add image of detector
        fig.add_trace(go.Heatmap(z=stack_image, colorscale=dcolorsc, x=x, y=y, hoverinfo='text',
                                 colorbar={'tickvals':tickvals, 'ticktext':ticktext,
                                           'title':{'text':'# times flagged as CR', 'side':'right', 'font':{'size':18}}},
                                 name=''))

        # add extraction box
        fig.add_trace(go.Scatter(x=np.arange(len(box_upper)), y=box_upper, mode="lines",
            line=go.scatter.Line(color='#222222', dash='dash'), showlegend=False,
            opacity=0.7, line_shape='hv', name='extraction box'))
        fig.add_trace(go.Scatter(x=np.arange(len(box_lower)), y=box_lower, mode="lines",
            line=go.scatter.Line(color='#222222', dash='dash'), showlegend=False,
            opacity=0.7, line_shape='hv', name='extraction box'))

        # y-axis zoom ranges
        zoom_options = [{'label':'Full Detector', 'yaxis_range':[0, stack_shape[0]]},
                        {'label':'Extraction Box', 'yaxis_range':[(min(box_lower)-20), (max(box_upper)+20)]}]

        # Add the toggle buttons
        button_options = [{'label':zoom_options[0]['label'], 'method':'relayout',
                           'args':[{'yaxis.range': zoom_options[0]['yaxis_range']}]},
                          {'label':zoom_options[1]['label'], 'method':'relayout',
                           'args':[{'yaxis.range': zoom_options[1]['yaxis_range']}]}]

        fig.update_layout(updatemenus=[{'type'       : 'dropdown',
                                        'direction'  : 'down',
                                        'buttons'    : button_options,
                                        'pad'        : {'r':0, 't':0},
                                        'showactive' : True,
                                        'x'          : 1.07,  # Position of the buttons; might require some more tweaking
                                        'xanchor'    : 'right',
                                        'y'          : 1.07,
                                        'yanchor'    : 'top'}])

        # Set the initial y-axis range (Full View)
        fig.update_yaxes(range=[0, stack_shape[0]])
        fig.update_xaxes(range=[0, stack_shape[1]])

        fig.update_layout(width=stack_shape[1] + 50,
                          height=int(stack_shape[1] * stack_shape[1] / stack_shape[0]) ) # adds space for colorbar to not squeeze the x-axis
        fig.update_layout(title={'text':title_text, 'x':0.5}, font={'family':'Arial, sans-serif', 'size':16})
        fig.write_html(file_path)


def split_plot(splits, box_lower, box_upper, split_num, individual_exposure_times, texpt,
               obs_id, propid, plot_dir, interactive):
    """Creates a visualization of where CR pixels are in each subexposure

    Parameters
    ----------
    splits: list
        list of CR flagged pixels in each subexposure

    box_lower: array
        1d array of ints of the bottom of the extraction box (0 indexed)

    box_upper: array
        1d array of ints of the top of the extraction box (0 indexed)

    split_num: int
        Number of splits in the stack

    individual_exposure_times: list
        List of exposure times for each subexposure

    texpt: float
        Value of total exposure time

    obs_id: str
        ipppssoot of observation

    propid: int
        proposal id of observation

    plot_dir: str
        Directory to save plot in

    interactive: bool 
        If True, uses Plotly to create an interactive zoomable html plot
    """
    custom_cmap = colors.ListedColormap([
        'k', 'tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:cyan', 'tab:olive',
        'tab:purple', 'tab:pink', 'tab:brown', 'tab:grey',])
    cmap = colors.ListedColormap(_gen_color(custom_cmap, 3))

    bounds = np.arange(4)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # Define grid, dependent on number of splits:
    if ((len(splits)) % 2) == 0:
        nrows = len(splits) / 2
    else:
        nrows = (len(splits) + 1) / 2

    row_value = int(nrows)

    if not interactive:
        fig, ax = plt.subplots(nrows=row_value, ncols=2, figsize=(9, nrows * 2))
        ax = ax.flatten()

        # Plot each subexposure with CR pixels a different color
        for num, axis in enumerate(ax):
            if num < len(splits):
                axis.imshow(splits[num], interpolation='none', origin='lower',
                    extent=(0, splits[num].shape[1], 0, splits[num].shape[0]),
                    cmap=cmap, norm=norm, aspect='auto')
                axis.step(np.arange(len(box_upper)), box_upper, color='#222222', where='post',
                    lw=0.7, alpha=0.7, ls='--')
                axis.step(np.arange(len(box_lower)), box_lower, color='#222222', where='post',
                    lw=0.7, alpha=0.7, ls='--')

                if ((splits[num].shape[0] - max(box_upper)) > 20) and (min(box_lower) > 20):
                    axis.set_ylim([min(box_lower) - 20, max(box_upper) + 20])
                    axis.set_title(f"zoomed subexposure {(num+1)!s}, exposure time {individual_exposure_times[num]!s}")

                else:
                    axis.set_title(f"subexposure {(num + 1)!s}, exposure time {individual_exposure_times[num]!s}")

            else:
                axis.set_axis_off()

        fig.suptitle(f"CR flagged pixels in individual splits for: {obs_id}\n Proposal {propid!s}, " \
                     f"total exposure time {texpt:.2f}, {split_num!s} subexposures")
        fig.tight_layout()

        plot_name = obs_id + '_splits.png'
        file_path = os.path.join(plot_dir, plot_name)
        plt.savefig(file_path, dpi=150, bbox_inches='tight')
        plt.close()

    else:
        subplot_titles = [f'zoomed subexposure {i+1}, exposure time {individual_exposure_times[i]}'
                          for i in range(len(splits))]

        title_text = f"CR flagged pixels in individual splits for: {obs_id}<br>Proposal {propid!s}, " \
                     f"total exposure time {texpt:.2f}, {split_num!s} subexposures"
        plot_name = f"{obs_id}_splits.html"
        file_path = os.path.join(plot_dir, plot_name)

        # Make Plotly figure
        fig = make_subplots(row_value, 2, horizontal_spacing=0.15, subplot_titles=subplot_titles)

        # Set up discrete color values
        dcolorsc = _discrete_colorscale(bvals=list(bounds[:-1]), colors=cmap.colors[:-1])

        # Add plots in each subplot
        row_iterator = 1
        for num, split in enumerate(splits):
            # calculate required x and y range to not center the pixels at 0,0
            x = np.arange(start=0, stop=split.shape[1]+1, step=1)
            y = np.arange(start=0, stop=split.shape[0]+1, step=1)

            # determine correct row, column to put the plot in
            if (num + 1) % 2 != 0:
                current_row = row_iterator
                row_iterator += 1
                current_column = 1
            else:
                #current_row = current_row
                current_column = 2

            # plot the pixel of each split and the extraction box values
            fig.add_trace(go.Heatmap(z=split, colorscale=dcolorsc, showscale=False, x=x, y=y, hoverinfo='text'),
                          current_row, current_column)
            fig.add_trace(go.Scatter(x=np.arange(len(box_upper)), y=box_upper, mode="lines",
                                     line=go.scatter.Line(color='#222222', dash='dash'),
                                     showlegend=False, opacity=0.7, line_shape='hv', name='extraction box'),
                          current_row, current_column)
            fig.add_trace(go.Scatter(x=np.arange(len(box_lower)), y=box_lower, mode="lines",
                                     line=go.scatter.Line(color='#222222', dash='dash'), showlegend=False,
                                     opacity=0.7, line_shape='hv', name='extraction box'),
                          current_row, current_column)

            # zoom the plot to near the extraction region
            fig.update_yaxes(range=[min(box_lower) - 20, max(box_upper) + 20])
            fig.update_xaxes(range=[0, split.shape[1]])

            # make plots zoom at the same time
            fig.update_xaxes(matches='x')
            fig.update_yaxes(matches='y')

        fig.update_layout(width=splits[-1].shape[1] * 1.75,
                          height=(((max(box_upper) + 20) - (min(box_lower) - 20)) * 3 * len(splits)))
        fig.update_layout(title={'text': title_text, 'x': 0.5, 'y': 1. - (0.15 / len(splits))},
                          font={'family': 'Arial, sans-serif', 'size': 16},
                          title_pad={'b': 20*len(splits)})
        fig.write_html(file_path)


def call_ocrreject_exam():
    """Command line usage of ocrreject_exam
    """
    parser = argparse.ArgumentParser(
        description="Calculate fractions of cosmic ray rejected pixels inside and outside " \
                    "of an extraction box to test for CR algorithm failures.",
        epilog=f'v{__version__};  Written by {__author__}')

    parser.add_argument(dest='obs_ids', metavar='obs_id', type=str, nargs='+',
        help='observation id(s) in ipppssoot format')
    parser.add_argument('-d', dest='data_dir', type=str, default='.',
        help="directory containing observation flt and sx1/x1d files. Defaults to current " \
             "working directory.")
    parser.add_argument('-p', dest='plot', action='store_true',
        help="option to create diagnostic plots")
    parser.add_argument('-o', dest='plot_dir', type=str, default=None,
        help="output directory to store diagnostic plots if plot=True. Defaults to data_dir.")
    parser.add_argument('-i', dest='interactive', action='store_true',
        help="option to create zoomable html plots instead of static pngs. Defaults to False "
             "and requires Plotly if True")

    kwargs = vars(parser.parse_args())
    kwargs['verbose']=True
    ocrreject_exam(**kwargs)


if __name__ == '__main__':
    call_ocrreject_exam()
