#! /usr/bin/env python
import os
import numpy as np
import astropy.io.fits as fits
import argparse

__author__ = 'Joleen K. Carlberg & Matt Dallas'
__version__ = 1.0

def ocrrej_exam(obsid, dir=None):
    ''' Compares the rate of cosmic rays in the extraction box and everywhere else 
        in a CCD spectroscopic image. Based on crrej_exam from STIS ISR 2019-02.

        Higher ratios of cosmic ray rates in the extraction box to the rest of the image
        may indicate the need to rerun stistools.ocrreject() with different parameters.

        Parameters
        ----------
        obsid : str
            A STIS observation ID rootname in ipppssoots format (ie odvkl1040)

        dir : str
            Directory containing both the flat fielded (_flt.fits) and extracted 
            spectrum (_sx1.fits or _x1d.fits) files of the observation.

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
            
        If called from the command line, prints the avg extraction, outside, and ratio values for quick verification 

    '''
    if not dir:
        dir = os.getcwd()+'/'

    # Get flt and sx1 filepaths
    flt_file = os.path.join(dir, obs_id+'_flt.fits')

    if not os.path.exists(flt_file):
      raise IOError(f"No _flt file in working directory for {obs_id}")

    sx1_file = os.path.join(dir, obs_id+'_sx1.fits')

    if not os.path.exists(sx1_file):
      raise IOError(f"No _flt file in working directory for {obs_id}")

    # Check that the number of sci extensions matches the number of crsplits
    with fits.open(flt_file) as flt_hdul:
        nrptexp_num = flt_hdul[0].header['NRPTEXP']
        crsplit_num = flt_hdul[0].header['CRSPLIT']
        sci_num = len([hdu.name for hdu in flt_hdul if "SCI" in hdu.name]) # Counts the number of sci extensions

    if ((crsplit_num)*(nrptexp_num))-(sci_num)!= 0:
        raise ValueError(f"cr-split or nrptexp value in flt header does not match the number of sci extentsions for {obs_id}")

    # Calculate cr fraction in and out of extraction box
    with fits.open(sx1_file) as sx1_hdul:
        spec = sx1_hdul[1].data[0]
        shdr = sx1_hdul[0].header

    extrlocy = spec['EXTRLOCY']-1 # y coord floats of the middle of the extraction box 
    del_pix = spec['EXTRSIZE']/2. # float value the extraction box extends above or below extrlocy
    box_upper=np.ceil(extrlocy+del_pix).astype(int) # Ints of pixel values above end of trace bc python is upper bound exclusive 
    box_lower=np.floor(extrlocy-del_pix).astype(int) # Ints of pixel values below end of trace

    # Fill each of these lists with values for each cr split
    extr_fracs = [] # float of the fraction of pixels flagged as cr inside the extraction box
    outside_fracs = [] # float of the fraction of pixels flagged as cr outside the extraction box
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

    return results

def call_ocrrej_exam():
    '''Command line usage of ocrrejectexam'''
    
    parser = argparse.ArgumentParser(description='Calculate fractions of cosmic ray rejected pixels inside and outside of an extraction box to test for cr algorithm failures.',
        epilog=f'v{__version__};  Written by {__author__}')

    parser.add_argument(dest='obsids', nargs='*', help='observation ids in ipppssoots format')
    parser.add_argument('-d', dest='dir', default=None, 
        help="directory containing observation flt and sx1 files os.getcwd()+'/'")
    args = parser.parse_args()
    
    for obsid in args.obsids:    
       print(f'Analyzing {osbid}:')
       result = ocrrej_exam(obsid, dir=args.dir)
    
       print('Fraction of Pixels Rejected as CRs')
       print(f"  Average across all extraction boxes: {result['avg_extr_frac']:.1%}")
       print(f"  Average across all external regions: {result['avg_outside_frac']:.1%}")
       print(f"  Average ratio between the two: {result['avg_ratio']:.1%}")


if __name__ == '__main__':
    call_crrej_exam()