import astropy
from astropy.io import fits
from astropy.modeling import models
import numpy as np
import matplotlib
from matplotlib import pyplot
import os
import sys

file1=sys.argv[1]

subcube=fits.open(file1)
data=subcube[0].data
header=subcube[0].header
print(data.shape)

rm_cube=data[:,:,:]

rm_cube[rm_cube==0]=np.nan

rms_cube=np.nanstd(rm_cube, axis=0)

#snr=rm_cube/rms_cube

#subcube[0].data=np.float32(snr)
#subcube.writeto(file1+"snr.fits", clobber=True)
subcube[0].data=np.float32(rms_cube)
subcube.writeto(file1+"_rms.fits", clobber=True)
