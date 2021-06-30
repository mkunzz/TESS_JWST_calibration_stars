#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 09:19:58 2021

@author: mkunz
"""
import lightkurve as lk
import numpy as np
from astroquery.mast import Catalogs
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------

radSearch = 4/60 #radius in degrees

catalogData = Catalogs.query_object('TIC 219820925', radius = radSearch, catalog = "TIC")
ra = catalogData[0]['ra']
dec = catalogData[0]['dec']

# Print out the first row in the table
catalogData[:5]['ID', 'Tmag', 'Jmag', 'ra', 'dec', 'objType']

coord = SkyCoord(ra, dec, unit = "deg")

sectorTable = Tesscut.get_sectors(coord)


#-----------------------------------------------------------------------------
# Define a function to simplify the plotting command that we do repeatedly.
def plot_cutout(image):
    #
    #Plot image and add grid lines.
    #
    plt.imshow(image, origin = 'lower', cmap = plt.cm.YlGnBu_r, 
           vmax = np.percentile(image, 92),
           vmin = np.percentile(image, 5))

    plt.grid(axis = 'both',color = 'white', ls = 'solid')

def aperture_phot(image, aperture):
    #
    #Sum-up the pixels that are in the aperture for one image.
    #image and aperture are 2D arrays that need to be the same size.
    #
    #aperture is a boolean array where True means to include the light of those pixels.
    
    flux = np.sum(image[aperture])

    return flux

def make_lc(flux_data, aperture):
    #
    #Apply the 2d aperture array to the and time series of 2D images. 
    #Return the photometric series by summing over the pixels that are in the aperture.
    #
    #Aperture is a boolean array where True means it is in the desired aperture.
    
    
    flux = np.array(list (map (lambda x: aperture_phot(x, aperture), flux_data) ) )

    return flux
#-----------------------------------------------------------------------------

hdulist = Tesscut.get_cutouts(coord, 12)

print('number of sectors = ' + str(len(hdulist))) #number of sectors
i = 0
while i < len(hdulist):
    print('\n')
    print(i)
    hdu1 = hdulist[i]

    # Use all pixels in our aperture.
    aperture = hdu1[2].data == 1
    flux1 = make_lc(hdu1[1].data['FLUX'], aperture)
    time1 = hdu1[1].data['TIME']

    # Plot the flux change of the dimmest pixels by using percentile.
    bkgAperture = hdu1[1].data['FLUX'][0] < np.percentile(hdu1[1].data['FLUX'][0], 5)
    bkgFlux1 = make_lc(hdu1[1].data['FLUX'], bkgAperture)
    
    # subtract background
    bkgSubFlux = flux1 - (bkgFlux1 * np.sum(aperture) / np.sum(bkgAperture) )
    normalized = bkgSubFlux/np.mean(bkgSubFlux)
    
    #plotting lightcurve of tesscut, normalized and background subtracted
    plt.scatter(time1, normalized, s=.8, c='k')
    y_lower_lim = np.mean(normalized)-1.5*np.std(normalized)
    y_upper_lim = np.mean(normalized)+1.5*np.std(normalized)
    plt.ylim(y_lower_lim, y_upper_lim)
    plt.xlabel('TIME (BTJD)')
    plt.ylabel('Flux (e-/s)')
    plt.title('Background Subtracted Flux, Normalized')
    plt.show()
    plt.close()
    #plt.savefig('/Users/mkunz/tesscuts/' + 'TIC 219820925_' + str(i) + 'lc' + '.png')

    i+=1


'''
tc = lk.search_tesscut('TIC 219820925').download_all()


lc = tc[0].to_lightcurve()
#lc2 = lc.flatten(window_length=1501)

flux = lc.flux.value
time = lc.time.value

bkg = flux < np.percentile(flux, 5)
bkg_subtracted = flux - bkg

#y_lower_lim = np.mean(bkg_subtracted)-np.std(bkg_subtracted)
#y_upper_lim = np.mean(bkg_subtracted)+np.std(bkg_subtracted)

plt.scatter(time, bkg_subtracted, s=.8, c='k')
#plt.ylim(y_lower_lim, y_upper_lim)
#plt.ylim(.99, 1.01)
plt.xlabel('time [days]')
plt.ylabel('flux [e- s-]')
plt.title('lightcurve of TIC 219820925')
'''