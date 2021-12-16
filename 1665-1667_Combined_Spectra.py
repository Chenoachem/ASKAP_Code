#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
#get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


# In[4]:
Cube_1665=sys.argv[1]
Cube_1667=sys.argv[2]
RA = sys.argv[3]
DEC = sys.argv[4]

datacube_1665 = fits.open(Cube_1665)
data_1665 = datacube_1665[0].data
header_1665 = datacube_1665[0].header
data_1665.shape


# In[6]:


datacube_1667 = fits.open(Cube_1667)
data_1667 = datacube_1667[0].data
header_1667 = datacube_1667[0].header
data_1667.shape


# In[7]:


rp = datacube_1665[0].header['CRPIX3']
rf = datacube_1665[0].header['CRVAL3']
df = datacube_1665[0].header['CDELT3']
nf = datacube_1665[0].header['NAXIS3']
xvals_1665 = rf + df*(np.arange(nf)-rp)
#xvals are the frequency in Hz np.subtract(xvals,1.66555e+09)
vels_1665=np.multiply(np.subtract(1.665402e+09,xvals_1665),0.180012068)
#Correct for m/s into km/s
vels_1665=np.divide(vels_1665,1000)
#print(vels_1665[1])


# In[8]:


rp_2 = datacube_1667[0].header['CRPIX3']
rf_2 = datacube_1667[0].header['CRVAL3']
df_2 = datacube_1667[0].header['CDELT3']
nf_2 = datacube_1667[0].header['NAXIS3']
xvals_1667 = rf_2 + df_2*(np.arange(nf_2)-rp_2)
vels_1667=np.multiply(np.subtract(1.667539e+09,xvals_1667),0.1797813772)
#Correct for m/s into km/s
vels_1667=np.divide(vels_1667,1000)
#print(vels_1667[1])


# In[28]:


ra=RA
dec=DEC
c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
r=c.ra.hms
d=c.dec.dms
GAL=c.galactic
c_gal=c.transform_to('galactic')

gala=c_gal.to_string('decimal', precision=3)
#print(gala)


# In[29]:


w_1665 = wcs.WCS(header_1665, naxis=2)
xpix_1665,ypix_1665=c.to_pixel(w_1665,origin=0,mode='wcs')
xpix_1665=int(xpix_1665)
ypix_1665=int(ypix_1665)

#print(xpix_1665)
#print(ypix_1665)


# In[30]:


w_1667 = wcs.WCS(header_1667, naxis=2)
xpix_1667,ypix_1667=c.to_pixel(w_1667,origin=0,mode='wcs')
xpix_1667=int(xpix_1667)
ypix_1667=int(ypix_1667)

#print(xpix_1667)
#print(ypix_1667)


# In[31]:


signal_1665=[]
for x in range(0, 2998):
    value_1665 = np.nanmean(data_1665[x,ypix_1665:ypix_1665+1,xpix_1665:xpix_1665+1])
    #print value
    signal_1665.append(value_1665)


# In[32]:


signal_1667=[]
for x in range(0, 2998):
    value_1667 = np.nanmean(data_1667[x,ypix_1667:ypix_1667+1,xpix_1667:xpix_1667+1])
    #print value
    signal_1667.append(value_1667)


# In[55]:


max_signal_1665=np.nan_to_num(signal_1665)
max_value_1665 = np.amax(max_signal_1665)
sliced_1665=np.argmax(max_signal_1665, axis=0)
max_signal_1667=np.nan_to_num(signal_1667)
max_value_1667 = np.amax(max_signal_1667)
sliced_1667=np.argmax(max_signal_1667, axis=0)
#print(sliced_1665,sliced_1667)
#a=sliced_1665 - abs(sliced_1665-sliced_1667)-200
#b=sliced_1665 + abs(sliced_1665-sliced_1667)+100
#sliced_1665=1416
a=sliced_1665-200
b=sliced_1665+200


# In[59]:


import matplotlib.patches as mpatches
#Set the minor axis counters
majorLocator = MultipleLocator(5)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(0.5)
        
#Make a spectra        
bigfig_1=plt.figure(figsize=(12,5))
ax1=bigfig_1.add_subplot(111)
ax1.step(vels_1665[0:2998],signal_1665[0:2998],color='black',label='1665')
ax1.step(vels_1665[0:2998],signal_1667[0:2998],color='red', label='1667')
#ax1.step(xvals[0:2000],signal[0:2000],color='red')
ax1.set_title("G"+gala, fontsize=24)
ax1.set_xlabel("Velocity (km/s)",fontsize=12)
ax1.set_xlim(vels_1665[a],vels_1665[b])
ax1.set_ylabel("Flux Density (Jy)",fontsize=12)
ax1.tick_params(labelsize=12, labelbottom=True)
ax1.ticklabel_format(useOffset=False)
ax1.xaxis.set_major_locator(majorLocator)
ax1.xaxis.set_major_formatter(majorFormatter)
#plum_patch = mpatches.Patch(color='grey')
#ax1.fill_between(vels[0:3000], mean_sig-RMS, mean_sig+RMS, facecolor='grey', alpha=0.5)
ax1.legend(fontsize=12)
plt.grid(True)

# for the minor ticks, use no labels; default NullFormatter
ax1.xaxis.set_minor_locator(minorLocator)

#Save the figure.     
bigfig_1.savefig("G"+gala+"_oneplot.png")


# In[60]:


import matplotlib.patches as mpatches
fig, axs = plt.subplots(2, sharex=True, sharey=False)
axs[0].step(vels_1665[0:2998],signal_1665[0:2998],color='black', label='1665')
axs[0].set_title("G"+gala, fontsize=12)
axs[0].set_xlim(vels_1665[a],vels_1665[b])
axs[0].set_ylabel("Flux Density (Jy)",fontsize=12)
axs[0].tick_params(labelsize=12, labelbottom=False)
axs[0].legend(fontsize=12)

axs[1].step(vels_1665[0:2998],signal_1667[0:2998],color='red',label='1667')
axs[1].set_xlabel("LSR Velocity (km/s)",fontsize=12)
axs[1].set_xlim(vels_1665[a],vels_1665[b])
axs[1].set_ylabel("Flux Density (Jy)",fontsize=12)
axs[1].tick_params(labelsize=12, labelbottom=True)
axs[1].legend(fontsize=12)

axs[1].xaxis.set_minor_locator(minorLocator)
plt.tight_layout()
#fig.savefig("Circular_Polarisation_RASP_Maser_Plot.pdf", dpi=300)
fig.savefig("G"+gala+"_twoplot_Samex.png")


# In[64]:


#import matplotlib.patches as mpatches
#fig, axs = plt.subplots(2, sharex=False, sharey=True)
#axs[0].step(vels_1665[0:2998],signal_1665[0:2998],color='black', label='1665')
#axs[0].set_title("G"+gala, fontsize=12)
#axs[0].set_xlim(vels_1665[sliced_1665-200],vels_1665[sliced_1665+200])
#axs[0].set_ylabel("Flux Density (Jy)",fontsize=12)
#axs[0].tick_params(labelsize=12, labelbottom=True)
#axs[0].legend(fontsize=12)

#axs[1].step(vels_1667[0:2998],signal_1667[0:2998],color='red',label='1667')
#axs[1].set_xlabel("LSR Velocity (km/s)",fontsize=12)
#axs[1].set_xlim(vels_1667[sliced_1667-200],vels_1667[sliced_1667+200])
#axs[1].set_ylabel("Flux Density (Jy)",fontsize=12)
#axs[1].tick_params(labelsize=12, labelbottom=True)
#axs[1].legend(fontsize=12)
#plt.tight_layout()
#axs[1].xaxis.set_minor_locator(minorLocator)
#fig.savefig("G"+gala+"_twoplot_Diffx.png")


# In[ ]:




