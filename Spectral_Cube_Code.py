#!/usr/bin/env python
import os
import argparse
import numpy as np
import astropy
from astropy.io import fits
import matplotlib
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy import wcs
import math
import sys
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from spectral_cube import SpectralCube

mycube='image.restored.i.SB23442.cube.notaper.rob_0.5.maxUV_0.GASKAP_GP_OH_A.beam01.fits'
out='GASKAP_OH_Beam_Sub.fits'

OH_data = fits.open(mycube)
cube = SpectralCube.read(OH_data)
header = OH_data[0].header

ra='16:35:34'
dec='-47:31:11'
c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
w = wcs.WCS(header, naxis=2)
xpix,ypix=c.to_pixel(w,origin=0,mode='wcs')
xpix=int(xpix)
ypix=int(ypix)

scube = cube[:,ypix-30:ypix+30,xpix-30:xpix+30]
scube.write(out, overwrite=True)
OH_data.close()

small_cube = SpectralCube.read(out)
OH_small=fits.open(out)
header=OH_small[0].header
ra='16:35:34'
dec='-47:31:11'

c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
w = wcs.WCS(header, naxis=2)
xpix,ypix=c.to_pixel(w,origin=0,mode='wcs')
xpix=int(xpix)
ypix=int(ypix)
print(ypix,xpix)

max_signal=np.amax(small_cube[:,ypix,xpix])
print(max_signal)


max = np.argmax(small_cube[:,ypix,xpix], axis=0)
max=int(max)
print max

#cube2 = small_cube.with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=1665.402 * u.MHz)

#print cube2

#sub_cube_slab = cube2.spectral_slab(-150. *u.km / u.s, -50. *u.km / u.s)

plt.plot(sub_cube_slab[:,ypix,xpix])


plt.savefig('spectra_test.png')


