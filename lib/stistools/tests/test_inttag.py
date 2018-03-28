from ..stistools.inttag import inttag
from astropy.io.fits import printdiff

def iraf_compare_results1():
    iraf_file = 'data/'
    #inttag(iraf_file, iraf_file.replace('.fits', '_out1.fits'))
    inttag()

def iraf_compare_results2():
    iraf_file = 'data/'
    inttag()

def iraf_compare_results3():
    iraf_file = 'data/'
    inttag()