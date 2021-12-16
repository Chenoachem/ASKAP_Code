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
from astropy.wcs import utils
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from spectral_cube import SpectralCube
from astropy.wcs.utils import pixel_to_skycoord
import aplpy
matplotlib.use('Agg')
from astropy.io.votable import parse_single_table


mycube=sys.argv[1]
out=sys.argv[2]


datacube = fits.open(mycube)
data = datacube[0].data
header = datacube[0].header
#print(data.shape)

#print datacube[0].header['CUNIT3']
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
#print(vels)

w = wcs.WCS(header, naxis=2)
ra1, dec1 = w.all_pix2world(15, 15, 1, ra_dec_order=True)
c=SkyCoord(ra1,dec1,unit='deg',frame='fk5')
r=c.ra.hms
d=c.dec.dms

signal=[]
for x in range(0, 3842):
    value = np.nanmedian(data[x,15:15+1,15:15+1])
    #print value
    signal.append(value)

max_signal=np.nan_to_num(signal)
max_value = np.amax(max_signal)
sliced=np.argmax(max_signal, axis=0)
a=sliced-100
b=sliced+100

#print(vels[sliced])

rms_number=np.nanstd(data[sliced,15+8:15+12,15+8:15+12])
sig1=rms_number*3
sig2=rms_number*6
sig3=rms_number*9
sig4=rms_number*12
sig5=rms_number*15
sig6=rms_number*18
#print(sig4)


file = open("Sources_2.txt", "a")
file.write("\n"+"RA ="+ str(ra1) + "\n"+"Dec ="+ str(dec1) + "\n"+"Peak Intensity = "+str(max_signal[sliced]) + "\n"+"Velocity ="+ str(vels[sliced]))
file.close

#Set the minor axis counters
majorLocator = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1)

sli=int(sliced)

        
#Make a spectra
bigfig=plt.figure(figsize=(11,5))

fig=aplpy.FITSFigure(mycube,figure=bigfig,dimensions=[0, 1],slices=[sli,],subplot=(1,2,1))
#fig.recenter(ra1,dec1,width=0.25,height=0.25)
#fig.show_contour(levels=(0.6,0.9,1.2,1.5,1.8), colors='black', linewidths=1)
fig.show_contour(levels=(sig1,sig2,sig3,sig4,sig5,sig6), colors='black', linewidths=1)
fig.add_beam()
fig.beam.show(zorder=100)
fig.beam.set_corner('bottom right')
fig.beam.set_frame(True)
fig.beam.set_hatch('/')
fig.beam.set_color('black')
fig.beam.set_edgecolor('black')
fig.beam.set_facecolor('black')

fig.tick_labels.set_xformat("hh:mm:ss")
fig.tick_labels.set_yformat("dd:mm:ss")
fig.set_title(("RA: "+str(np.rint(r[0]))+":"+str(np.rint(r[1]))+":"+str(np.round(r[2],2))+" " + "Dec"+str(np.rint(d[0]))+":"+str(np.rint(abs(d[1])))+":"+str(np.round(abs(d[2])))), fontsize='14')
#fig.show_markers(posloc[0], posloc[1], edgecolor='red', linewidth=2, marker='*', s=s_plot*(fov_wise/fov_mwa)**2)
fig.set_tick_color('black')


ax1=bigfig.add_subplot(122)
ax1.step(vels,signal,color='black')
#ax1.set_title("OH Emission Source"+str(np.round(ra1,2)) +","+str(np.round(dec1,2)) , fontsize=12)
ax1.set_xlabel("Velocity (km/s)",fontsize=12)
ax1.set_xlim(vels[a],vels[b])
ax1.set_ylabel("Intensity (Jy/beam)",fontsize=12)
ax1.tick_params(labelsize=12, labelbottom=True)
ax1.ticklabel_format(useOffset=False)
#ax1.xaxis.set_major_locator(majorLocator)
#ax1.xaxis.set_major_formatter(majorFormatter)

# for the minor ticks, use no labels; default NullFormatter
ax1.xaxis.set_minor_locator(minorLocator)

plt.tight_layout(pad=4, w_pad=4, h_pad=4)

#Save the figure.
bigfig.savefig(out)
