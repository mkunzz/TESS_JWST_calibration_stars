#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 13:21:18 2021

@author: mkunz
"""
import lightkurve as lk
import numpy as np
from astroquery.mast import Catalogs
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
#from astropy.io import fits
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------
# Search for JWST Calibration stars using
# lightkurve for Tess 2 min data using TIC numbers

candles_with_TIC = ['TIC 327587572', 'TIC 247923021', 'TIC 149505899', 'TIC 352817378', 'TIC 147921014', 'TIC 471015233', 'TIC 383553764', 'TIC 166698220', 'TIC 41232189', 'TIC 298165335', 'TIC 198456033', 'TIC 219820925', 'TIC 181240911', 'TIC 80313923', 'TIC 75586606', 'TIC 165370459', 'TIC 229945862', 'TIC 440765193', 'TIC 219752116', 'TIC 135656809', 'TIC 27533327', 'TIC 39464221', 'TIC 441120034', 'TIC 397558558', 'TIC 140282069', 'TIC 420814525', 'TIC 32869782', 'TIC 54837036', 'TIC 365653206', 'TIC 229980646', 'TIC 8591766', 'TIC 417544924', 'TIC 144599609', 'TIC 315627636', 'TIC 207440438', 'TIC 219094190', 'TIC 233095291', 'TIC 219114641', 'TIC 233067231', 'TIC 233075513', 'TIC 219897252', 'TIC 233205654']

candles_without_2min = []
i = 0
while i < len(candles_with_TIC):
    lc = lk.search_lightcurve(candles_with_TIC[i], author='SPOC', exptime=120)
    if len(lc) == 0:
        candles_without_2min.append(candles_with_TIC[i])
    else:
        pass
    i+=1

print('Length of candles with TIC = ' + str(len(candles_with_TIC)))
print('Length of candles without 2 minute data = ' + str(len(candles_without_2min)))
print('-----------------------------------------------------------------------')
#-----------------------------------------------------------------------------


#print(candles_without_2min)

#-----------------------------------------------------------------------------

#none = ['TIC 352817378']
#has_FFI =['TIC 219820925']
none = []
has_FFI =[]
RAs = []
DECs = []
coordinates = []
sectorTables = []
nearby_Stars =[]
i = 0
while i < len(candles_without_2min):
    starName = candles_without_2min[i]
    radSearch = 4/60 #radius in degrees
    
    catalogData = Catalogs.query_object(starName, radius = radSearch, catalog = "TIC")
    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']
    RAs.append(ra)
    DECs.append(dec)
    # Print out the first row in the table
    #print( catalogData[:5]['ID', 'Tmag', 'Jmag', 'ra', 'dec', 'objType'] )

    coord = SkyCoord(ra, dec, unit = "deg")
    coordinates.append(coord)
    sectorTable = Tesscut.get_sectors(coord)
    sectorTables.append(sectorTable)
    # Create a list of nearby bright stars (tess magnitude less than 14) from the rest of the data for later.
    bright = catalogData['Tmag'] < 15
    
    # Make it a list of Ra, Dec pairs of the bright ones. This is now a list of nearby bright stars.
    nearbyStars = list( map( lambda x,y:[x,y], catalogData[bright]['ra'], catalogData[bright]['dec'] ) )
    len(nearbyStars)
    nearby_Stars.append(nearbyStars)
    
    if len(sectorTable) == 0:
        none.append(candles_without_2min[i])
    else:
        has_FFI.append(candles_without_2min[i])
        
    i+=1


#print(none)
#print(has_FFI)
#-----------------------------------------------------------------------------
# Define a function to simplify the plotting command that we do repeatedly.
def plot_cutout(image):
    """
    Plot image and add grid lines.
    """
    plt.imshow(image, origin = 'lower', cmap = plt.cm.YlGnBu_r, 
           vmax = np.percentile(image, 92),
           vmin = np.percentile(image, 5))

    plt.grid(axis = 'both',color = 'white', ls = 'solid')

def aperture_phot(image, aperture):
    """
    Sum-up the pixels that are in the aperture for one image.
    image and aperture are 2D arrays that need to be the same size.
    
    aperture is a boolean array where True means to include the light of those pixels.
    """
    flux = np.sum(image[aperture])

    return flux

def make_lc(flux_data, aperture):
    """
    Apply the 2d aperture array to the and time series of 2D images. 
    Return the photometric series by summing over the pixels that are in the aperture.
    
    Aperture is a boolean array where True means it is in the desired aperture.
    """
    
    flux = np.array(list (map (lambda x: aperture_phot(x, aperture), flux_data) ) )

    return flux
#-----------------------------------------------------------------------------
i = 0
while i < len(has_FFI):
    
    hdulist = Tesscut.get_cutouts(coordinates[i], 20)
    
    hdulist[0].info()
    hdulist[0][0].header['SECTOR']

    hdu1 = hdulist[0]
    firstImage = hdu1[1].data['FLUX'][0]
    
    fig = plt.figure(figsize=(7, 7))
    plot_cutout(firstImage)
    plt.xlabel('Image Column',fontsize = 14)
    plt.ylabel('Image Row',fontsize = 14)
    
    
    hdu2 = hdulist[1]

    firstImage = hdu2[1].data['FLUX'][0]
    
    wcs = WCS(hdu2[2].header)
    
    fig = plt.figure(figsize = (8, 8))
    fig.add_subplot(111, projection = wcs)
    plot_cutout(firstImage)
    
    plt.xlabel('RA', fontsize = 12)
    plt.ylabel('Dec', fontsize = 12)
    
    
    starloc = wcs.all_world2pix([[RAs[i],DECs[i]]],0)  #Second is origin
    plt.scatter(starloc[0,0], starloc[0,1],s = 45,color = 'red')
    
    # Plot nearby stars as well, which we created using our Catalog call above.
    nearbyLoc = wcs.all_world2pix(nearby_Stars[i][1:],0)
    plt.scatter(nearbyLoc[1:, 0], nearbyLoc[1:, 1], 
                s = 25, color = 'orange')
    
    # Use all pixels in our aperture.
    aperture = hdu1[2].data == 1
    
    flux1 = make_lc(hdu1[1].data['FLUX'], aperture)
    time1 = hdu1[1].data['TIME']
    
    plt.figure(figsize = (11,5))
    plt.plot(time1, flux1, 'k.-', lw = .5)
    plt.xlim(1325,1342)
    plt.xlabel('TIME (BTJD)')
    plt.ylabel('Flux (e-/s)')
    plt.title('Flux in Photometric Aperture')
    
    # Plot the flux change of the dimmest pixels by using percentile.
    bkgAperture = hdu1[1].data['FLUX'][0] < np.percentile(hdu1[1].data['FLUX'][0], 5)
    
    bkgFlux1 = make_lc(hdu1[1].data['FLUX'], bkgAperture)
    time1 = hdu1[1].data['TIME']
    
    plt.figure(figsize = (11, 5))
    plt.plot(time1, bkgFlux1, 'r.-', lw = .5)
    
    plt.xlim(1325, 1342)
    plt.xlabel('TIME (BTJD)')
    plt.ylabel('Estimate of Background')
    plt.title('Background Flux')
    
    bkgSubFlux = flux1 - (bkgFlux1 * np.sum(aperture) / np.sum(bkgAperture) )

    plt.figure(figsize = (11,5))
    plt.plot(time1, bkgSubFlux,'.-k', lw = 0.5)
    
    plt.xlim(1325, 1336)
    plt.xlabel('TIME (BTJD)')
    plt.ylabel('Flux (e-/s)')
    plt.title('Background Subtracted Flux')
    
    i+=1








































