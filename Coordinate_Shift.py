#!/usr/bin/env python
import os
import argparse
import numpy as np
import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy import wcs
import math
import sys
from spectral_cube import SpectralCube

mycube=sys.argv[1]
#out='GASKAP_OH_Beam_Sub.fits'
cat_files=sys.argv[2]

#OH_data = fits.open(mycube)
#cube = SpectralCube.read(OH_data)
#header = OH_data[0].header

file=open(cat_files)
for line in file:
    RA,DEC = line.strip().split(",")
    RA=str(RA)
    DEC=str(DEC)
    #c = SkyCoord(ra=RA, dec=DEC, unit='deg', frame='fk5')
    
#ra='16:35:34'
#dec='-47:31:11'
    c = SkyCoord(ra=RA, dec=DEC, unit=(u.hourangle, u.deg))
    GAL=c.galactic
    c_gal=c.transform_to('galactic')
    gala=c_gal.to_string('decimal')
    w = wcs.WCS(header, naxis=2)
    #xpix,ypix=c.to_pixel(w,origin=0,mode='wcs')
    #xpix=int(xpix)
    #ypix=int(ypix)
    #max_signal=np.amax(cube[:,ypix,xpix])
    #max_signal=str(max_signal)
    

    #scube = cube[:,ypix-15:ypix+15,xpix-15:xpix+15]
    #scube.write("Cubelet"+"_"+str(xpix)+"_"+str(ypix)+".fits", overwrite=True)
    
    
    file = open("Sources_1667.txt", "a")
    file.write("\n" + "RA = " + RA + "," +"Dec = "+DEC + ","+"Gal = "+gala + ","  +" )
    file.close

OH_data.close()
    #file = open("Sources.txt", "a")
    #file.write(ra,dec,max_signal+ "\n")
    #file.close()



