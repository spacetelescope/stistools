#! /usr/bin/env python
import os
import warnings
import numpy as np
import astropy.io.fits as fits
import argparse
import matplotlib.cm as colormap
import matplotlib.colors as colors
import matplotlib.pyplot as plt

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

USER_WARNED = False

__author__ = 'Joleen K. Carlberg & Matt Dallas'
__version__ = 1.0

def ocrreject_exam(obs_id=None, data_dir=None, flt=None, sx1=None, plot=False, plot_dir=None, interactive=False):
    """Compares the rate of cosmic rays in the extraction box and everywhere else 
    in a CCD spectroscopic image. Based on crrej_exam from STIS ISR 2019-02.

    Higher ratios of cosmic ray rates in the extraction box to the rest of the image
    may indicate the need to rerun stistools.ocrreject() with different parameters.

    Parameters
    ----------
    obsid : str
        A STIS observation ID rootname in ipppssoot format (ie odvkl1040)

    data_dir : str
        Directory containing both the flat fielded (_flt.fits) and extracted 
        spectrum (_sx1.fits or _x1d.fits) files of the observation. 
        Defaults to pwd and requires trailing /

    flt : str
        Path to flt file. Useful if flt and sx1 are in different locations or have custom names.

    sx1 : str
        Path to sx1 file. Useful if flt and sx1 are in different locations or have custom names.

    plot : bool
        Option to generate diagnostic plots, default=False

    plot_dir : str
        Directory to save diagnostic plots in if plot=True. 
        Defaults to dir parameter and requires trailing /

    interactive : bool
        Option to generate zoomable html plots using plotly, default=False

    Returns
    -------
    results : dict
        Dictionary containing 
            extr_fracs : cr rejection rates in the extraction boxes for each crsplit
            outside_fracs : cr rejection rates outside the extraction boxes for each crsplit
            ratios : extr_fracs/outside_fracs
            avg_extr_frac : The average of extr_fracs
            avg_outside_frac : The average of outside_fracs
            avg_ratio : avg_extr_frac/avg_outside_frac
        
    If called from the command line, prints the avg extraction, outside, and ratio values for quick verification.

    """
    global USER_WARNED

    if not data_dir:
        data_dir = os.getcwd()+'/'

    if not plot_dir:
        plot_dir = data_dir

    if obs_id is None:
        if flt is None or sx1 is None:
            raise ValueError("If 'obs_id' is not provided, both 'flt' and 'sx1' must be specified.")
        else:
            flt_file = flt
            sx1_file = sx1

    elif obs_id is not None:
        if flt is not None or sx1 is not None:
            raise ValueError("If 'obs_id' is provided, both 'flt' and 'sx1' must not be provided.")
        else:
            # Get flt and sx1/x1d filepaths
            flt_file = os.path.join(data_dir, obs_id+'_flt.fits')

            if not os.path.exists(flt_file):
                raise IOError(f"No _flt file in {data_dir} for {obs_id}")

            sx1_file = os.path.join(data_dir, obs_id+'_sx1.fits')

            if not os.path.exists(sx1_file):
                sx1_file = os.path.join(data_dir, obs_id+'_x1d.fits') # if sx1 doesn't exist check for custom made x1d
                if not os.path.exists(sx1_file):
                    raise IOError(f"No _sx1 file in {data_dir} for {obs_id}")

    # Check that the number of sci extensions matches the number of crsplits
    with fits.open(flt_file) as flt_hdul:
        propid = flt_hdul[0].header['PROPOSID']
        rootname = flt_hdul[0].header['ROOTNAME']
        nrptexp_num = flt_hdul[0].header['NRPTEXP']
        crsplit_num = flt_hdul[0].header['CRSPLIT']
        sci_num = len([hdu.name for hdu in flt_hdul if "SCI" in hdu.name]) # Counts the number of sci extensions

    if ((crsplit_num)*(nrptexp_num))-(sci_num)!= 0:
        raise ValueError(f"cr-split or nrptexp value in flt header does not match the number of sci extentsions for {obs_id}")

    # Calculate cr fraction in and out of extraction box
    with fits.open(sx1_file) as sx1_hdul:
        spec = sx1_hdul[1].data[0]
        shdr = sx1_hdul[0].header

    extrlocy = spec['EXTRLOCY']-1 # y coords of the middle of the extraction box 
    del_pix = spec['EXTRSIZE']/2. # value the extraction box extends above or below extrlocy
    box_upper=np.ceil(extrlocy+del_pix).astype(int) # Ints of pixel values above end of trace bc python is upper bound exclusive 
    box_lower=np.floor(extrlocy-del_pix).astype(int) # Ints of pixel values below end of trace

    # Fill each of these lists with values for each cr split
    extr_fracs = [] # fraction of pixels flagged as cr inside the extraction box for each split
    outside_fracs = [] # fraction of pixels flagged as cr outside the extraction box for each split
    cr_rejected_locs = [] # 2d array of 1s where a cr exists and 0 elsewhere
    exposure_times = []

    with fits.open(flt_file) as flt_hdul:
        flt_shape = flt_hdul['sci', 1].data.shape # shape of the data 
        
        # Check that the extraction box doesn't extend beyond the image: this breaks the method
        if np.any(box_lower < 0) or np.any(box_upper-1 > flt_shape[0]): # Subtract 1 because the box extends to the value of the pixel before
            raise ValueError(f"Extraction box coords extend above or below the cr subexposures for {obs_id}")

        extr_mask = np.zeros(flt_shape)
        outside_mask = np.ones(flt_shape)

        for column in range(0,flt_shape[1]):
            extr_mask[box_lower[column]:box_upper[column], column] = 1 # 1s inside the extraction box, 0s outside
            outside_mask[box_lower[column]:box_upper[column], column] = 0 # 0s inside the extraction box, 1s outside

        n_extr = np.count_nonzero(extr_mask) # number of pixels inside the extraction box
        n_outside = np.count_nonzero(outside_mask) # number of pixels outside the extraction box

        for i, hdu in enumerate(flt_hdul):
            if hdu.name == 'SCI':
                exposure_times.append(hdu.header['EXPTIME'])
                dq_array = flt_hdul[i+2].data # dq array corresponding to each sci extentsion

                extr_rej_pix = np.zeros(flt_shape) # 2d array where there is a 1 if a pixel inside the extraction box is marked as a cr
                np.place(extr_rej_pix, ((extr_mask == 1) & (dq_array & 2**13 != 0)), 1)

                outside_rej_pix = np.zeros(flt_shape) # 2d array where there is a 1 if a pixel outside the extraction box is marked as a cr
                np.place(outside_rej_pix, ((outside_mask == 1) & (dq_array & 2**13 != 0)), 1)

                extr_cr_count = np.count_nonzero(extr_rej_pix)
                outside_cr_count = np.count_nonzero(outside_rej_pix)

                extr_fracs.append(extr_cr_count/n_extr)
                outside_fracs.append(outside_cr_count/n_outside)

                cr_rejected_pix = extr_rej_pix+outside_rej_pix
                cr_rejected_locs.append(cr_rejected_pix)

    extr_fracs = np.asarray(extr_fracs)
    outside_fracs = np.asarray(outside_fracs)
    ratios = extr_fracs/outside_fracs # ratio of extraction to outside the box in each image

    avg_extr_frac = (np.sum(extr_fracs))/(len(extr_fracs)) # Average fraction of crs inside extraction box
    avg_outside_frac = (np.sum(outside_fracs))/(len(outside_fracs)) # Average fraction of crs outside extraction box
    avg_ratio = avg_extr_frac/avg_outside_frac # Average ratio of the stack

    results ={'extr_fracs':extr_fracs, 'outside_fracs':outside_fracs, 'ratios':ratios, 'avg_extr_frac':avg_extr_frac, 'avg_outside_frac':avg_outside_frac, 'avg_ratio':avg_ratio}

    if plot and not interactive: # case with interactive = false
        cr_rejected_stack = np.sum(cr_rejected_locs, axis=0) # stack all located crs on top of eachother
        stacked_exposure_time = sum(exposure_times)
        stack_plot(cr_rejected_stack, box_lower, box_upper, len(cr_rejected_locs), stacked_exposure_time, rootname, propid, plot_dir, interactive=interactive)
        split_plot(cr_rejected_locs, box_lower, box_upper, len(cr_rejected_locs), exposure_times, stacked_exposure_time, rootname, propid, plot_dir, interactive=interactive)
    
    elif plot and interactive and HAS_PLOTLY: # case with interactive = True and plotly is installed
        cr_rejected_stack = np.sum(cr_rejected_locs, axis=0) # stack all located crs on top of eachother
        stacked_exposure_time = sum(exposure_times)
        stack_plot(cr_rejected_stack, box_lower, box_upper, len(cr_rejected_locs), stacked_exposure_time, rootname, propid, plot_dir, interactive=interactive)
        split_plot(cr_rejected_locs, box_lower, box_upper, len(cr_rejected_locs), exposure_times, stacked_exposure_time, rootname, propid, plot_dir, interactive=interactive)
    
    elif plot and interactive and not USER_WARNED: # case with interactive = True and plotly is not installed
        warnings.warn('Plotly required for intercative plotting')
        USER_WARNED = True

    return results

# Plotting specific functions:
def gen_color(cmap, n):
    """Generates n distinct colors from a given colormap. 
    
    Based on mycolorpy's gen_color() from https://github.com/binodbhttr/mycolorpy"""

    colorlist = []

    for c in cmap.colors[0:n]:
        clr = colors.rgb2hex(c) # convert to hex
        colorlist.append(str(clr)) # create a list of these colors
    
    colorlist.pop(0) # Make it light grey rather than black at the beginning (I think it's easier on the eyes)
    colorlist.insert(0, '#F5F5F5')
        
    return colorlist

def discrete_colorscale(bvals, colors):
    """Takes desired boundary values and colors from a matplotlib colorplot and makes a plotly colorscale.
    
    Based on discrete_colorscale() from https://community.plot.ly/t/colors-for-discrete-ranges-in-heatmaps/7780"""

    if len(bvals) != len(colors)+1:
        raise ValueError('len(boundary values) should be equal to  len(colors)+1')
    bvals = sorted(bvals)     
    nvals = [(v-bvals[0])/(bvals[-1]-bvals[0]) for v in bvals]  # normalized values
    
    dcolorscale = [] # discrete colorscale
    for k in range(len(colors)):
        dcolorscale.extend([[nvals[k], colors[k]], [nvals[k+1], colors[k]]])
        
    return dcolorscale    

def generate_intervals(n, divisions):
    """Creates a list of strings that are the positions requred for centering an evenly spaced colorbar in plotly"""

    result = np.linspace(0, n, divisions, endpoint=False)
    offset = (result[1]-result[0])/2
    result = result + offset
    result = [str(x) for x in list(result)[:len(list(result))]]
    
    return result

def stack_plot(stack_image, box_lower, box_upper, split_num, texpt, obs_id, propid, plot_dir, interactive):
    """Creates a visualization of where cr pixels are in a stacked image  

    Parameters
    ----------
    stack_image : array
        2d array to plot.

    box_lower : array
        1d array of ints of the bottom of the extraction box 0 indexed. 

    box_upper : array
        1d array of ints of the top of the extraction box 0 indexed.

    split_num : int
        Number of splits in the stack.

    texpt : float
        Value of total exposure time.

    obs_id : str
        ipppssoot of observation

    propid : int
        proposal id of observation

    plot_dir : str
        Directory to save plot in. Requires trailing /

    interactive : bool 
        If True, uses plotly to create an interactive zoomable html plot
    """
    
    stack_shape = stack_image.shape
    max_stack_value = int(np.max(stack_image)) # This is usually equal to stack_shape,
    # in the case where a cr pixel is not in all splits at the same location this value should be used
    custom_cmap = colors.ListedColormap(['k', 'tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:cyan', 'tab:olive', 'tab:purple', 'tab:pink', 'tab:brown', 'tab:grey'])
    cmap = colors.ListedColormap(gen_color(custom_cmap, max_stack_value+1))
    bounds = np.arange(max_stack_value+2)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    if not interactive:
        # create matplotlib image
        fig, (ax1,ax2,ax3) = plt.subplots(nrows=1, ncols=3, figsize=(9,20*(9/41)), gridspec_kw={'width_ratios': [1, 1, 0.05], 'height_ratios': [1]})

        for axis in [ax1,ax2]:
            axis.imshow(stack_image, interpolation='none', origin="lower", extent=(0, stack_shape[1], 0, stack_shape[0]), cmap=cmap, norm=norm, aspect='auto')
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

        cb = fig.colorbar(colormap.ScalarMappable(norm=norm, cmap=cmap), cax=ax3, label='# times flagged as cr', ticks=np.arange(max_stack_value, max_stack_value+2)-0.5)
        cb.set_ticklabels(np.arange(max_stack_value, max_stack_value+2)-1)

        fig.suptitle('CR flagged pixels in stacked image: '+obs_id+'\n Proposal '+str(propid)+', exposure time '+f'{texpt:.2f}'+', '+str(split_num)+' subexposures')
        fig.tight_layout()

        plot_name = obs_id + '_stacked.png'
        file_path = plot_dir + plot_name
        plt.savefig(file_path, dpi=150, bbox_inches='tight')
        plt.close()
    
    else:
        # Create plotly image
        fig = go.Figure()

        # calculate required x and y range, colorbar info, and figure titles
        x = np.arange(start=0, stop=stack_shape[1]+1, step=1)
        y = np.arange(start=0, stop=stack_shape[0]+1, step=1)

        dcolorsc = discrete_colorscale(bvals=list(bounds), colors=cmap.colors)

        ticktext = [str(x) for x in list(bounds)[:len(list(bounds))-1]]
        tickvals = generate_intervals(len(ticktext)-1, len(ticktext))

        title_text = 'CR flagged pixels in stacked image: '+obs_id+'<br>'+'Proposal '+str(propid)+', exposure time '+f'{texpt:.2f}'+', '+str(split_num)+' subexposures'
        plot_name = obs_id + '_stacked.html'
        file_path = plot_dir + plot_name

        # add image of detector
        fig.add_trace(go.Heatmap(z=stack_image, colorscale=dcolorsc, x=x, y=y, hoverinfo='text', colorbar={'tickvals':tickvals, 'ticktext':ticktext, 'title':{'text':'# times flagged as cr', 'side':'right', 'font':{'size':18}}}, name=''))

        # add extraction box
        fig.add_trace(go.Scatter(x=np.arange(len(box_upper)),y=box_upper,mode="lines",line=go.scatter.Line(color='#222222', dash='dash'),showlegend=False, opacity=0.7, line_shape='hv', name='extraction box'))
        fig.add_trace(go.Scatter(x=np.arange(len(box_lower)),y=box_lower,mode="lines",line=go.scatter.Line(color='#222222', dash='dash'),showlegend=False, opacity=0.7, line_shape='hv',  name='extraction box'))

        # y-axis zoom ranges
        zoom_options = [{'label':'Full Detector', 'yaxis_range':[0, stack_shape[0]]}, 
                        {'label':'Extraction Box', 'yaxis_range':[(min(box_lower)-20), (max(box_upper)+20)]}]

        # Add the toggle buttons
        button_options = [{'label':zoom_options[0]['label'], 'method':'relayout', 'args':[{'yaxis.range':zoom_options[0]['yaxis_range']}]},
                        {'label':zoom_options[1]['label'], 'method':'relayout', 'args':[{'yaxis.range': zoom_options[1]['yaxis_range']}]}]

        fig.update_layout(updatemenus=[{'type':'dropdown',
                                        'direction':'down',
                                        'buttons':button_options,
                                        'pad':{'r':0, 't':0},
                                        'showactive':True,
                                        'x': 1.07,  # Position of the buttons- also might require some more tweaking
                                        'xanchor':'right',
                                        'y': 1.07,
                                        'yanchor':'top'}])

        # Set the initial y-axis range (Full View)
        fig.update_yaxes(range=[0, stack_shape[0]])
        fig.update_xaxes(range=[0, stack_shape[1]])

        fig.update_layout(width=stack_shape[1]+ 50, height=int(stack_shape[1] * stack_shape[1] / stack_shape[0]) ) # adds space for colorbar to not squeeze the x axis
        fig.update_layout(title={'text':title_text, 'x':0.5}, font={'family':'Arial, sans-serif', 'size':16})
        fig.write_html(file_path)

def split_plot(splits, box_lower, box_upper, split_num, individual_exposure_times, texpt, obs_id, propid, plot_dir, interactive):
    """Creates a visualization of where cr pixels are in each subexposure  

    Parameters
    ----------
    splits : list
        list of cr placements in each subexposure (ie the cr_rejected_locs output of ocrreject_exam)

    box_lower : array
        1d array of ints of the bottom of the extraction box 0 indexed. 

    box_upper : array
        1d array of ints of the top of the extraction box 0 indexed.

    split_num : int
        Number of splits in the stack, (ie len(cr_rejected_locs)).

    individual_exposure_times: list
        List of exposure times for each subexposure

    texpt : float
        Value of total exposure time

    obs_id : str
        ipppssoot of observation

    propid : int
        proposal id of observation

    plot_dir : str
        Directory to save plot in. Requires trailing /

    interactive : bool 
        If True, uses plotly to create an interactive zoomable html plot
    """

    custom_cmap = colors.ListedColormap(['k', 'tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:cyan', 'tab:olive', 'tab:purple', 'tab:pink', 'tab:brown', 'tab:grey'])
    cmap = colors.ListedColormap(gen_color(custom_cmap, 3))

    bounds = np.arange(4)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # Define grid, dependent on number of splits:
    if ((len(splits))%2) == 0:
        nrows =  (len(splits))/2
    else:
        nrows = ((len(splits))+1)/2

    row_value = int(nrows)

    if not interactive:
        fig, ax = plt.subplots(nrows=row_value, ncols=2, figsize=(9, nrows*2))
        ax = ax.flatten()

        # Plot each subexposure with cr pixels a different color
        for num, axis in enumerate(ax):
            if num<len(splits):
                axis.imshow(splits[num], interpolation='none', origin='lower', 
                            extent=(0, splits[num].shape[1], 0, splits[num].shape[0]),
                            cmap=cmap, norm=norm, aspect='auto')
                axis.step(np.arange(len(box_upper)), box_upper, color='#222222', where='post', lw=0.7, alpha=0.7, ls='--')
                axis.step(np.arange(len(box_lower)), box_lower, color='#222222', where='post', lw=0.7, alpha=0.7, ls='--')

                if ((splits[num].shape[0] - max(box_upper)) > 20) and (min(box_lower) >20):
                    axis.set_ylim([(min(box_lower)-20),(max(box_upper)+20)])
                    axis.set_title('zoomed subexposure '+str(num+1)+', exposure time '+str(individual_exposure_times[num]))

                else:
                    axis.set_title('subexposure '+str(num+1)+', exposure time '+str(individual_exposure_times[num]))

            else:
                axis.set_axis_off()

        fig.suptitle('CR flagged pixels in individual splits for: '+obs_id+ '\n Proposal '+str(propid)+', total exposure time '+f'{texpt:.2f}'+', '+str(split_num)+' subexposures')
        fig.tight_layout()

        plot_name = obs_id + '_splits.png'
        file_path = plot_dir + plot_name
        plt.savefig(file_path, dpi=150, bbox_inches='tight')
        plt.close()
    
    else:
        subplot_titles = [f'zoomed subexposure {i+1}, exposure time {individual_exposure_times[i]}' for i in range(len(splits))]

        title_text = 'CR flagged pixels in individual splits for: '+obs_id+ '<br>'+'Proposal '+str(propid)+', total exposure time '+f'{texpt:.2f}'+', '+str(split_num)+' subexposures'
        plot_name = obs_id + '_splits.html'
        file_path = plot_dir + plot_name

        # Make plotly figure
        fig = make_subplots(row_value, 2, horizontal_spacing=0.15, subplot_titles=subplot_titles)

        # Set up discrete color values
        dcolorsc = discrete_colorscale(bvals=list(bounds[:-1]), colors=cmap.colors[:-1])

        # Add plots in each subplot
        row_iterator = 1
        for num, split in enumerate(splits):
            # calculate required x and y range to not center the pixels at 0,0
            x = np.arange(start=0, stop=split.shape[1]+1, step=1)
            y = np.arange(start=0, stop=split.shape[0]+1, step=1)
            
            # determine correct row, column to put the plot in
            if (num+1)%2 != 0:
                current_row = row_iterator
                row_iterator+=1
            else:
                current_row = current_row
            
            if (num+1)%2 == 0:
                current_column = 2
            else:
                current_column = 1

            # plot the pixel of each split and the extraction box values
            fig.add_trace(go.Heatmap(z=split, colorscale=dcolorsc, showscale=False, x=x, y=y, hoverinfo='text'), current_row, current_column)
            fig.add_trace(go.Scatter(x=np.arange(len(box_upper)),y=box_upper,mode="lines",line=go.scatter.Line(color='#222222', dash='dash'),showlegend=False, opacity=0.7, line_shape='hv', name='extraction box'), current_row, current_column)
            fig.add_trace(go.Scatter(x=np.arange(len(box_lower)),y=box_lower,mode="lines",line=go.scatter.Line(color='#222222', dash='dash'),showlegend=False, opacity=0.7, line_shape='hv',  name='extraction box'), current_row, current_column)

            # zoom the plot to near the extraction region
            fig.update_yaxes(range=[min(box_lower)-20,max(box_upper)+20])
            fig.update_xaxes(range=[0, split.shape[1]])

            # make plots zoom at the same time
            fig.update_xaxes(matches='x')
            fig.update_yaxes(matches='y')
        
        fig.update_layout(width=split.shape[1]*1.75, height=((((max(box_upper)+20)- (min(box_lower)-20))*3*len(splits))))
        fig.update_layout(title={'text':title_text, 'x':0.5, 'y':1.0-(0.15/len(splits))}, font={'family':'Arial, sans-serif', 'size':16}, title_pad={'b': 20*len(splits)})
        fig.write_html(file_path)

def call_ocrreject_exam():
    """Command line usage of ocrreject_exam"""
    
    parser = argparse.ArgumentParser(description='Calculate fractions of cosmic ray rejected pixels inside and outside of an extraction box to test for cr algorithm failures.',
        epilog=f'v{__version__};  Written by {__author__}')

    parser.add_argument('--obs', dest='obs_ids', nargs='*', default=None, help='observation ids in ipppssoots format')
    parser.add_argument('--flt', dest='flt', default=None, help='path to flt file')
    parser.add_argument('--sx1', dest='sx1', default=None, help='path to sx1 file')
    parser.add_argument('--d', dest='data_dir', default=None, help="directory containing observation flt and sx1 files. Defaults to pwd and requires trailing /")
    parser.add_argument('-p', dest='plot', help="option to create diagnostic plots", action='store_true')
    parser.add_argument('--pd', dest='plot_dir', default=None, help="directory to store diagnostic plots if plot=True. Defaults to data_dir argument and requires trailing /")
    parser.add_argument('-i', dest='interactive', default=False, help="option to create zoomable html plots instead of static pngs. Defaults to False and requires plotly if True")

    args = parser.parse_args()
    
    if args.obs_ids is not None:
        if args.flt is not None or args.sx1 is not None:
            raise ValueError("If 'obs_id' is provided, both 'flt' and 'sx1' must not be provided.")
        else:
            for obsid in args.obs_ids:    
                print(f'\nAnalyzing {obsid}:')
                result = ocrreject_exam(obsid=args.obs_ids, data_dir=args.data_dir, flt=args.flt, sx1=args.sx1, plot=args.plot, plot_dir=args.plot_dir, interactive=args.interactive)
                
                print('Fraction of Pixels Rejected as CRs')
                print(f"  Average across all extraction boxes: {result['avg_extr_frac']:.1%}")
                print(f"  Average across all external regions: {result['avg_outside_frac']:.1%}")
                print(f"  Average ratio between the two: {result['avg_ratio']:.2f}")

    elif args.flt is None or args.flt is None:
        raise ValueError("If 'obs_id' is not provided, both 'flt' and 'sx1' must be specified.")

    else:
        result = ocrreject_exam(obsid=args.obs_ids, data_dir=args.data_dir, flt=args.flt, sx1=args.sx1, plot=args.plot, plot_dir=args.plot_dir, interactive=args.interactive)
            
        print('Fraction of Pixels Rejected as CRs')
        print(f"  Average across all extraction boxes: {result['avg_extr_frac']:.1%}")
        print(f"  Average across all external regions: {result['avg_outside_frac']:.1%}")
        print(f"  Average ratio between the two: {result['avg_ratio']:.2f}")

    


if __name__ == '__main__':
    call_ocrreject_exam()