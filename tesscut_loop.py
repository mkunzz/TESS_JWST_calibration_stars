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
radSearch = 4/60 #radius in degrees

none = []
has_FFI =[]
i = 0
while i < len(candles_without_2min):
    starName = candles_without_2min[i]
    
    catalogData = Catalogs.query_object(starName, radius = radSearch, catalog = "TIC")
    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']
    
    coord = SkyCoord(ra, dec, unit = "deg")
    sectorTable = Tesscut.get_sectors(coord)
    
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
has_FFI = ['TIC 219820925',
 'TIC 80313923',
 'TIC 229945862',
 'TIC 440765193',
 'TIC 8591766',
 'TIC 417544924',
 'TIC 144599609',
 'TIC 207440438',
 'TIC 219094190',
 'TIC 233095291',
 'TIC 219114641',
 'TIC 233067231',
 'TIC 233075513',
 'TIC 219897252',
 'TIC 233205654']

    
catalogData = Catalogs.query_object(has_FFI[0], radius = radSearch, catalog = "TIC")
ra = catalogData[0]['ra']
dec = catalogData[0]['dec']

# Print out the first row in the table
catalogData[:5]['ID', 'Tmag', 'Jmag', 'ra', 'dec', 'objType']

coord = SkyCoord(ra, dec, unit = "deg")

sectorTable = Tesscut.get_sectors(coord)


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
    bkgAperture = hdu1[1].data['FLUX'][0] < np.nanpercentile(hdu1[1].data['FLUX'][0], 15)
    bkgFlux1 = make_lc(hdu1[1].data['FLUX'], bkgAperture)
    
    # subtract background
    bkgSubFlux = flux1 - (bkgFlux1 * np.nansum(aperture) / np.nansum(bkgAperture) )
    normalized = bkgSubFlux/np.nanmean(bkgSubFlux)
    
    #plotting lightcurve of tesscut, normalized and background subtracted
    plt.scatter(time1, normalized, s=.8, c='k')
    y_lower_lim = np.nanmean(normalized)-1.5*np.nanstd(normalized)
    y_upper_lim = np.nanmean(normalized)+1.5*np.nanstd(normalized)
    plt.ylim(y_lower_lim, y_upper_lim)
    plt.xlabel('TIME (BTJD)')
    plt.ylabel('Flux (e-/s)')
    plt.title('Background Subtracted Flux, Normalized')
    plt.show()
    plt.close()
    #plt.savefig('/Users/mkunz/tesscuts/' + 'TIC 219820925_' + str(i) + 'lc' + '.png')




'''
ii = 0
while ii < len(has_FFI):
    
    catalogData = Catalogs.query_object(has_FFI[ii], radius = radSearch, catalog = "TIC")
    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']
    
    # Print out the first row in the table
    catalogData[:5]['ID', 'Tmag', 'Jmag', 'ra', 'dec', 'objType']
    
    coord = SkyCoord(ra, dec, unit = "deg")
    
    sectorTable = Tesscut.get_sectors(coord)
    
    
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
        bkgAperture = hdu1[1].data['FLUX'][0] < np.nanpercentile(hdu1[1].data['FLUX'][0], 15)
        bkgFlux1 = make_lc(hdu1[1].data['FLUX'], bkgAperture)
        
        # subtract background
        bkgSubFlux = flux1 - (bkgFlux1 * np.nansum(aperture) / np.nansum(bkgAperture) )
        normalized = bkgSubFlux/np.nanmean(bkgSubFlux)
        
        #plotting lightcurve of tesscut, normalized and background subtracted
        plt.scatter(time1, normalized, s=.8, c='k')
        y_lower_lim = np.nanmean(normalized)-1.5*np.nanstd(normalized)
        y_upper_lim = np.nanmean(normalized)+1.5*np.nanstd(normalized)
        plt.ylim(y_lower_lim, y_upper_lim)
        plt.xlabel('TIME (BTJD)')
        plt.ylabel('Flux (e-/s)')
        plt.title('Background Subtracted Flux, Normalized')
        plt.show()
        plt.close()
        #plt.savefig('/Users/mkunz/tesscuts/' + 'TIC 219820925_' + str(i) + 'lc' + '.png')
    
        i+=1
    
    
    ii+=1
'''




































