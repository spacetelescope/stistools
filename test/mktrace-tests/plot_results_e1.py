"""
Plot results of test with E1 aperture positions 
input file: o8pp31020_crj.fits

"""

#!/usr/bin/env python

import pylab
import pyfits

###################
fref=pyfits.open('ref_o8pp31020_crj_1dt.fits)
a2disp_ref=fref[1].data.field('A2DISPL')
fref.close()

f=pyfits.open('o8pp31020_crj_1dt.fits)
a2disp=f[1].data.field('A2DISPL')
f.close()

ax1=pylab.subplot(231)
pylab.plot(a2disp[1792])
pylab.plot(a2disp_ref[1792])
pylab.title('trace tables')
###########################

fref=pyfits.open('ref_o8pp31020_crj_1dt__interp.fits)
data_ref=fref[0].data
fref.close()

f=pyfits.open('o8pp31020_crj_1dt_interp.fits)
data_ref=f[0].data
f.close()
ax1=pylab.subplot(232)
pylab.plot(data_ref)
pylab.plot(data)
pylab.title('interpolated trace')
##############################

fref=pyfits.open('ref_o8pp31020_crj_1dt__interpfit.fits)
data_ref=fref[0].data
fref.close()

f=pyfits.open('o8pp31020_crj_1dt_interpfit.fits)
data_ref=f[0].data
f.close()
ax1=pylab.subplot(233)
pylab.plot(data_ref)
pylab.plot(data)
pylab.title('fit to interp trace')


##############################

fref=pyfits.open('ref_o8pp31020_crj_1dt__sci.fits)
data_ref=fref[0].data
fref.close()

f=pyfits.open('o8pp31020_crj_1dt_sci.fits)
data_ref=f[0].data
f.close()
ax1=pylab.subplot(234)
pylab.plot(data_ref)
pylab.plot(data)
pylab.title('science trace')

##############################

fref=pyfits.open('ref_o8pp31020_crj_1dt__scifit.fits)
data_ref=fref[0].data
fref.close()

f=pyfits.open('o8pp31020_crj_1dt_scifit.fits)
data_ref=f[0].data
f.close()
ax1=pylab.subplot(235)
pylab.plot(data_ref)
pylab.plot(data)
pylab.title('fit to sci trace ')
pylab.show()
