#! /usr/bin/env python
import os
import numpy as np
import astropy.io.fits as fits
import argparse
from matplotlib import colormaps
import matplotlib.cm as colormap
import matplotlib.colors as colors
import matplotlib.pyplot as plt

__author__ = 'Joleen K. Carlberg & Matt Dallas'
__version__ = 1.0

def ocrreject_exam(obs_id, dir=None, plot=False, plot_dir=None):
    """Compares the rate of cosmic rays in the extraction box and everywhere else 
    in a CCD spectroscopic image. Based on crrej_exam from STIS ISR 2019-02.

    Higher ratios of cosmic ray rates in the extraction box to the rest of the image
    may indicate the need to rerun stistools.ocrreject() with different parameters.

    Parameters
    ----------
    obsid : str
        A STIS observation ID rootname in ipppssoot format (ie odvkl1040)

    dir : str
        Directory containing both the flat fielded (_flt.fits) and extracted 
        spectrum (_sx1.fits or _x1d.fits) files of the observation. 
        Defaults to pwd and requires trailing /

    plot : bool
        Option to generate diagnostic plots, default=False

    plot_dir : str
        Directory to save diagnostic plots in if plot=True. 
        Defaults to dir parameter and requires trailing /

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

    if not dir:
        dir = os.getcwd()+'/'

    if not plot_dir:
        plot_dir = dir

    # Get flt and sx1/x1d filepaths
    flt_file = os.path.join(dir, obs_id+'_flt.fits')

    if not os.path.exists(flt_file):
        raise IOError(f"No _flt file in {dir} for {obs_id}")

    sx1_file = os.path.join(dir, obs_id+'_sx1.fits')

    if not os.path.exists(sx1_file):
        sx1_file = os.path.join(dir, obs_id+'_x1d.fits') # if sx1 doesn't exist check for custom made x1d
        if not os.path.exists(sx1_file):
            raise IOError(f"No _sx1 file in {dir} for {obs_id}")

    # Check that the number of sci extensions matches the number of crsplits
    with fits.open(flt_file) as flt_hdul:
        propid = flt_hdul[0].header['PROPOSID']
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

    if plot:
        cr_rejected_stack = np.sum(cr_rejected_locs, axis=0) # stack all located crs on top of eachother
        stacked_exposure_time = sum(exposure_times)
        stack_plot(cr_rejected_stack, box_lower, box_upper, len(cr_rejected_locs), stacked_exposure_time, obs_id, propid, plot_dir)
        split_plot(cr_rejected_locs, box_lower, box_upper, len(cr_rejected_locs), exposure_times, stacked_exposure_time, obs_id, propid, plot_dir)

    return results

def stack_plot(stack_image, box_lower, box_upper, split_num, texpt, obs_id, propid, plot_dir):
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

    """

    stack_shape = stack_image.shape
    cmap = colors.ListedColormap(gen_color('turbo', split_num+1))
    bounds = np.arange(split_num+2)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, (ax1,ax2,ax3) = plt.subplots(nrows=1, ncols=3, figsize=(9,20*(9/41)), gridspec_kw={'width_ratios': [1, 1, 0.05], 'height_ratios': [1]})

    for axis in [ax1,ax2]:
        axis.imshow(stack_image, interpolation='nearest', origin="lower", extent=(0, stack_shape[1], 0, stack_shape[0]), cmap=cmap, norm=norm, aspect='auto')
        axis.step(np.arange(len(box_upper)), box_upper, color='w', where='post', lw=0.5, alpha=0.5, ls='--')
        axis.step(np.arange(len(box_lower)), box_lower, color='w', where='post', lw=0.5, alpha=0.5, ls='--')

    ax1.set_title('Full image')

    # If it is a large enough image, zoom the 2nd subplot around the extraction box region
    if ((stack_shape[0] - max(box_upper)) > 20) and (min(box_lower) > 20):
        ax2.set_ylim([(min(box_lower)-20),(max(box_upper)+20)])
        ax2.set_title('zoomed to 20 pixels above/below extraction box')

    # Otherwise just don't zoom in at all
    else:
        ax2.set_title('full image already 20 pixels above/below extraction box')

    fig.colorbar(colormap.ScalarMappable(norm=norm, cmap=cmap), cax=ax3, label='# times flagged as cr')

    fig.suptitle('CR flagged pixels in stacked image: '+obs_id+'\n Proposal '+str(propid)+', exposure time '+f'{texpt:.2f}'+', '+str(split_num)+' subexposures')
    fig.tight_layout()

    plot_name = obs_id + '_stacked.png'
    file_path = plot_dir + plot_name
    plt.savefig(file_path, dpi=150)
    plt.close()

def split_plot(splits, box_lower, box_upper, split_num, individual_exposure_times, texpt, obs_id, propid, plot_dir):
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
    """

    # Define grid, dependent on number of splits:
    if ((len(splits))%2) == 0:
        nrows =  (len(splits))/2
    else:
        nrows = ((len(splits))+1)/2

    row_value = int(nrows)

    fig, ax = plt.subplots(nrows=row_value, ncols=2, figsize=(9, nrows*2))
    ax = ax.flatten()

    cmap = colors.ListedColormap(gen_color('autumn', 3))
    bounds = np.arange(4)
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # Plot each subexposure with cr pixels a different color
    for num, axis in enumerate(ax):
        if num<len(splits):
            axis.imshow(splits[num], interpolation='nearest', origin='lower', 
                        extent=(0, splits[num].shape[1], 0, splits[num].shape[0]),
                        cmap=cmap, norm=norm, aspect='auto')
            axis.step(np.arange(len(box_upper)), box_upper, color='w', where='post', lw=0.5, alpha=0.5, ls='--')
            axis.step(np.arange(len(box_lower)), box_lower, color='w', where='post', lw=0.5, alpha=0.5, ls='--')

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
    plt.savefig(file_path, dpi=150)
    plt.close()

def gen_color(cmap, n):
    """Generates n distinct colors from a given colormap. 
    
    Based on mycolorpy's gen_color() from https://github.com/binodbhttr/mycolorpy"""

    c_map = colormaps[cmap]
    colorlist = []
    
    for c in np.linspace(0,1,n):
        rgba=c_map(c) # select the rgba value of the cmap at point c which is a number between 0 to 1
        clr=colors.rgb2hex(rgba) # convert to hex
        colorlist.append(str(clr)) # create a list of these colors

    colorlist.pop(0) # Make it dark grey rather than black at the beginning (I think it's easier on the eyes)
    colorlist.insert(0, '#A9A9A9')
        
    return colorlist

def call_ocrreject_exam():
    """Command line usage of ocrreject_exam"""
    
    parser = argparse.ArgumentParser(description='Calculate fractions of cosmic ray rejected pixels inside and outside of an extraction box to test for cr algorithm failures.',
        epilog=f'v{__version__};  Written by {__author__}')

    parser.add_argument(dest='obsids', nargs='*', help='observation ids in ipppssoots format')
    parser.add_argument('-d', dest='dir', default=None, help="directory containing observation flt and sx1 files. Defaults to pwd and requires trailing /")
    parser.add_argument('-p', dest='plot', help="option to create diagnostic plots", action='store_true')
    parser.add_argument('-pd', dest='plot_dir', default=None, help="directory to store diagnostic plots if plot=True. Defaults to dir argument and requires trailing /")
    
    args = parser.parse_args()
    #TODO add check that obsid arg actually has stuff stored in it
    
    for obsid in args.obsids:    
       print(f'\nAnalyzing {obsid}:')
       result = ocrreject_exam(obsid, dir=args.dir, plot=args.plot, plot_dir=args.plot_dir)
    
       print('Fraction of Pixels Rejected as CRs')
       print(f"  Average across all extraction boxes: {result['avg_extr_frac']:.1%}")
       print(f"  Average across all external regions: {result['avg_outside_frac']:.1%}")
       print(f"  Average ratio between the two: {result['avg_ratio']:.2f}")


if __name__ == '__main__':
    call_ocrreject_exam()