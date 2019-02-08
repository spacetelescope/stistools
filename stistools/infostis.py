#! /usr/bin/env python

import sys

__author__ = 'S. Tony Sohn'

def infostis(filenames):
    import numpy as np
    from astropy import units as u
    from astropy.coordinates import Angle
    from astropy.io import fits

    '''Print a header information for a list of STIS FITS images.

    Parameters:
        filenames : str
            Input STIS image name or list of STIS images.

    Output:
        Screen output with general hedaer infromation for STIS FITS images
        on a standard windows size of 80x24.
    '''

    if type(filenames) != list:
            filenames = [filenames]

    for file in filenames:
       with fits.open(file) as hdu:
          rootname = hdu[0].header['ROOTNAME'] # 
          proposid = hdu[0].header['PROPOSID'] # Proposal ID
          linenum  = hdu[0].header['LINENUM']  # Proposal logsheet line number = Exposure ID
          targname = hdu[0].header['TARGNAME'] #
          detector = hdu[0].header['DETECTOR'] #
          obstype  = hdu[0].header['OBSTYPE']  #
          obsmode  = hdu[0].header['OBSMODE']  #
          sclamp   = hdu[0].header['SCLAMP']   # lamp status, NONE or name of lamp which is on
          aperture = hdu[0].header['APERTURE'] #
          filtr    = hdu[0].header['FILTER']   #
          opt_elem = hdu[0].header['OPT_ELEM'] #

          ra_targ  = hdu[0].header['RA_TARG']  #
          dec_targ = hdu[0].header['DEC_TARG'] #
          equinox  = hdu[0].header['EQUINOX']  #
          binaxis1 = hdu[0].header['BINAXIS1'] # axis1 data bin size in unbinned detector pixels
          binaxis2 = hdu[0].header['BINAXIS2'] # axis2 data bin size in unbinned detector pixels
          texptime = hdu[0].header['TEXPTIME'] #
          nextend  = hdu[0].header['NEXTEND']  # Number of image sets combined during CR rejection

          subarray = hdu[0].header['SUBARRAY'] # T or F
          
          if obstype == 'SPECTROSCOPIC':
                  cenwave = hdu[0].header['CENWAVE']

          if detector == 'CCD':
                  ccdamp   = hdu[0].header['CCDAMP']   #
                  ccdgain  = hdu[0].header['CCDGAIN']  # 
                  crsplit  = hdu[0].header['CRSPLIT']  # 
    
          # Number of imsets
          if nextend > 0:
                  imsets = nextend / 3
          else:
                  imsets = 0

          rr = Angle(ra_targ, u.degree)
          rt = rr.hms
          dd = Angle(dec_targ, u.degree)
          dt = dd.signed_dms

          # Now the pretty printing part
          print("")
          print("-"*80)
          print("                                   S T I S")
          print("-"*80)
          print("")
          print("")
          print("%18s %-18s %16s %-25s" % ("Rootname:",rootname.upper(),"Detector:",detector))
          print("%18s %-18s %16s %-25s" % ("Proposal ID:",proposid,"Obs Type:",obstype))
          print("%18s %-18s %16s %-25s" % ("Exposure ID:",linenum,"Obs Mode:", obsmode))
          print("%54s %s"               % ("Lamp:",sclamp))
          print("%18s %-18s %16s %-25s" % ("Target Name:",targname,"Aperture:", aperture))
                    
          print("%19s %02.0f%s%02.0f%s%04.1f %23s %s" % ("Right Ascension: ", rt[0], ":", rt[1], ":", rt[2], "Filter:", filtr))
          
          if dt[0] >= 0:
                  print("%18s %s%02.0f%s%02.0f%s%04.1f %23s %s" % ("Declination:", "+", dt[1], ":", dt[2], ":", dt[3], "Opt Element:", opt_elem))
          else:
                  print("%18s %s%02.0f%s%02.0f%s%04.1f %23s %s" % ("Declination:", "-", dt[1], ":", dt[2], ":", dt[3], "Opt Element:", opt_elem))

          if obstype == 'SPECTROSCOPIC':
                  print("%18s %-18.1f %16s %-25d", "Equinox:" % (equinox,"Central Wave:", cenwave))
          elif detector == 'CCD':
                  print("%18s %-18.1f %16s %s" % ("Equinox:", equinox, "CCD amp:", ccdamp))
          else:
                  print("%18s %-18.1f" % ("Equinox:", equinox))

          if detector == 'CCD' and obstype == 'SPECTROSCOPIC':
                  print("%54s %s" % ("CCD amp:", ccdamp))
          elif detector == 'CCD':
                  print("%54s %d" % ("Gain:", ccdgain))

          if detector == 'CCD' and obstype == 'SPECTROSCOPIC':
                  print("%18s %-18d %16s %-25d" % ("Axis 1 binning:", binaxis1, "Gain:", ccdgain))
          elif detector == 'CCD':
                  print("%18s %-18d %16s %-25d" % ("Axis 1 binning:", binaxis1, "CR-split:", crsplit))
          else:
                  print("%18s %-18d" % ("Axis 1 binning:", binaxis1))

          if detector == 'CCD' and obstype == 'SPECTROSCOPIC':
                  print("%18s %-18d %16s %-25d" % ("Axis 2 binning:", binaxis2, "CR-split:", crsplit))
          else:
                  print("%18s %-18d" % ("Axis 2 binning:", binaxis2))
                  
          if subarray == 'T':
                  print("%18s %s" % ("Subarray:", "yes"))
          else:
                  print("%18s %s" % ("Subarray:", "no"))

          print("")
          print("%18s %-.1f %s" % ("Total Exp. Time:", texptime, "sec"))
          print("%18s %-18d" % ("Number of imsets:", imsets))
          print("")

    return      

def call_infostis(args):
    import argparse

    parser = argparse.ArgumentParser(
        description='Print a header information for a list of STIS FITS images.'),
    parser.add_argument(dest='filenames', metavar='filename', nargs='*', help='')
    args = parser.parse_args()

    infostis(args.filenames)

if __name__ == '__main__':
        call_infostis(sys.argv[1:])
        sys.exit(0)
