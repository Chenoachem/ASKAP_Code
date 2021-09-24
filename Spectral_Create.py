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

#In the terminal type in python Spectral_Create.py NameOfFitsCubelet NameofOutputSpectra

mycube=sys.argv[1]
out=sys.argv[2]

#Call in fits data cubelet
datacube = fits.open(mycube)
data = datacube[0].data
header = datacube[0].header
#print(data.shape)

#Read in and checked the frequency units
print datacube[0].header['CUNIT3']
rp = datacube[0].header['CRPIX3']
rf = datacube[0].header['CRVAL3']
df = datacube[0].header['CDELT3']
nf = datacube[0].header['NAXIS3']
xvals = rf + df*(np.arange(nf)-rp)
#xvals are the frequency in Hz np.subtract(xvals,1.66555e+09)
#xvals=xvals[300:999]
vels=np.multiply(np.subtract(1.665402e+09,xvals),0.180012)
#Correct for m/s into km/s
vels=np.divide(vels,1000)

#print(vels[100])


#ra='16:35:34'
#dec='-47:31:11'
#c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
#w = wcs.WCS(header, naxis=2)
#xpix,ypix=c.to_pixel(w,origin=0,mode='wcs')
#xpix=int(xpix)
#ypix=int(ypix)

#Create a numpy array with the flux information for each channel at pixel 30.  This value is set based on the centre of the make cubelet script.
signal=[]
for x in range(0, 3842):
    value = np.nanmean(data[x,30:30+1,30:30+1])
    #print value
    signal.append(value)

#Find the slice that contains the maximum signal.  Could be changed to min for absorption.
max_signal=np.nan_to_num(signal)
max_value = np.amax(max_signal)
sliced=np.argmax(max_signal, axis=0)
a=sliced-200
b=sliced+200

#print(vels[sliced])

#Open or create a text file named Sources where the inforamtion will be dropped to create a catalogue.
file = open("Sources.txt", "a")
file.write("\n"+"Peak Intensity = "+str(max_signal[sliced]) + "\n"+"Velocity ="+ str(vels[sliced]))
file.close

#Set the minor axis counters
majorLocator = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)
        
#Make a spectra
bigfig=plt.figure(figsize=(20,12))
ax1=bigfig.add_subplot(111)
ax1.step(vels,signal,color='blue')
ax1.set_title("OH Emission Source", fontsize=24)
ax1.set_xlabel("Velocity (km/s)",fontsize=24)
ax1.set_xlim(vels[a],vels[b])
ax1.set_ylabel("Intensity (Jy/beam)",fontsize=24)
ax1.tick_params(labelsize=22, labelbottom=True)
ax1.ticklabel_format(useOffset=False)
#ax1.xaxis.set_major_locator(majorLocator)
#ax1.xaxis.set_major_formatter(majorFormatter)

# for the minor ticks, use no labels; default NullFormatter
ax1.xaxis.set_minor_locator(minorLocator)

#Save the figure.
bigfig.savefig(out)
