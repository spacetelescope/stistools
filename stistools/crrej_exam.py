#! /usr/bin/env python
import glob
import numpy as np
import astropy.io.fits as fits
import os

__author__ = 'Joleen K. Carlberg'

def crrej_frac(obs_id, dir='./'):
   ''' Calculates the fraction of pixels that are cosmic ray rejected
       in the extraction region compared to the overal frame.
       
       Significantly higher rejection rates in the extraction region
       may require the user to tune the CR rejection parameters by hand.
       
       Parameters
       ----------
       obs_id : str
          A STIS observation ID, eg., odvlk1040
          
        dir : str (optional, default = ./)
          Directory containing both the flat field (_flt.fits) and extracted 
          spectrum (_sx1.fits or _x1d.fits) file of the observation.
          
        Returns
        -------
        frac_rej : float
          Fraction of pixels in the extracted region that are cosmic-ray rejected.

        frac_rej_expec : float
          The expected fraction of cosmic-ray rejected pixels (fraction computed for
           the entire CCD).
           
        Prints the rejection statistics to standard output.
   
   '''
   
   with fits.open(os.path.join(dir,obs_id+'_sx1.fits')) as spec_hdu:
         spec = spec_hdu[1].data[0]
         shdr =spec_hdu[0].header

   split_num = shdr['CRSPLIT']

   #Create a mask that defines where the extraction region is
   extr_mask = np.zeros((1024,1024))
   del_pix = spec["EXTRSIZE"]/2.
   for column in range(0,1024):
      row_mid =  spec['EXTRLOCY'][column]
      gd_row_low = int(np.ceil( row_mid - del_pix))
      gd_row_high = int(np.floor(row_mid + del_pix))
      extr_mask[gd_row_low:gd_row_high+1,column] = 1
      
   n_tot = np.count_nonzero(extr_mask) * split_num 
   
   with fits.open(os.path.join(dir,obs_id+'_flt.fits')) as flt_hdu:
        if len(flt_hdu) != 1+split_num*3:
           print('Unexpected number of extensions')
           
        n_rej = []
        
        for i in range(3,len(flt_hdu)+1 ,3):
           dat = flt_hdu[i].data
#           rej = dat[ (extr_mask == 1)  & (dat >= 8192)] #Data quality flag 8192 used for CR rejected pixels
           rej = dat[ (extr_mask == 1) & (dat & 2**13 != 0)] #Data quality flag 8192 (2^13) used for CR rejected pixels
           n_rej.append(np.count_nonzero(rej))
                 
   t_exp = float(shdr['TEXPTIME'])
   CR_rate = float(shdr['REJ_RATE'])
   
   # Calculate the rejection fraction and the rate of rejected pixels per sec
   n_pix_rej = np.sum(np.array(n_rej))
   frac_rej = n_pix_rej/float(n_tot)
   rej_rate = n_pix_rej * 1024.0**2 / (t_exp * n_tot)
   
   # Calculate the expected rejection fraction rate given the rate in header
   # This is the CR_Rate * ratio of number of pixels in full CCD vs extract region (n_tot /split_num)
   frac_rej_expec =  CR_rate * t_exp /(1024.*1024.*split_num)
   
   return frac_rej,frac_rej_expec

def call_crrej_exam():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Calculate fraction of cosmic rays flagged in extraction region')#, 
 #       epilog='{} package; Written by {}; v{}'.format(__package__, __author__, __version__))
    parser.add_argument(dest='rootnames', nargs='*', help='raw observations')
    parser.add_argument('-d', dest='dir', default='./', 
        help='data directory [default is current directory]')
    args = parser.parse_args()
    
    for root in args.rootnames:    
       print('Analyzing {}:'.format(root))
       rej_frac,expec_frac = crrej_frac(root, dir=args.dir)
    
       print('Percentage of Pixels Rejected as CRs')
       print('   Extraction Region: {:.1%}'.format(rej_frac))
       print('   Full CCD Frame: {:.1%} \n'.format(expec_frac))


if __name__ == '__main__':
    call_crrej_exam()

