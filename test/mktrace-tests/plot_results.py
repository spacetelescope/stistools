"""
Plot results of basic test level 1
input file: o45a03010_flt_crj.fits

"""

#!/usr/bin/env python

import pylab
import pyfits

###################
fref=pyfits.open('ref_o45a03010_flt_crj_1dt.fits')
a2disp_ref=fref[1].data.field('A2DISPL')
fref.close()

f=pyfits.open('o45a03010_flt_crj_1dt.fits')
a2disp=f[1].data.field('A2DISPL')
f.close()

pylab.subplot(231)
pylab.plot(a2disp[16])
pylab.plot(a2disp_ref[16], color='r')
pylab.title('trace tables')
###########################

fref=pyfits.open('ref_o45a03010_flt_crj_1dt_interp.fits')
data_ref=fref[0].data
fref.close()

f=pyfits.open('o45a03010_flt_crj_1dt_interp.fits')
data=f[0].data
f.close()
pylab.subplot(232)
pylab.plot(data_ref)
pylab.plot(data, color='r')
pylab.title('interpolated trace')
##############################

fref=pyfits.open('ref_o45a03010_flt_crj_1dt_interpfit.fits')
data_ref=fref[0].data
fref.close()

f=pyfits.open('o45a03010_flt_crj_1dt_interpfit.fits')
data=f[0].data
f.close()

pylab.subplot(233)
pylab.plot(data_ref)
pylab.plot(data, color='r')
pylab.title('fit to interp trace')


##############################

fref=pyfits.open('ref_o45a03010_flt_crj_1dt_sci.fits')
data_ref=fref[0].data
fref.close()

f=pyfits.open('o45a03010_flt_crj_1dt_sci.fits')
data=f[0].data
f.close()

pylab.subplot(234)
pylab.plot(data_ref)
pylab.plot(data, color='r')
pylab.title('science trace')

##############################

fref=pyfits.open('ref_o45a03010_flt_crj_1dt_scifit.fits')
data_ref=fref[0].data
fref.close()

f=pyfits.open('o45a03010_flt_crj_1dt_scifit.fits')
data=f[0].data
f.close()

pylab.subplot(235)
pylab.plot(data_ref)
pylab.plot(data,color='r')
pylab.title('fit to sci trace ')
pylab.show()
