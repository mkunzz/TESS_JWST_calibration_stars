#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 17:17:55 2021

@author: mkunz
"""

from astroquery.mast import Catalogs
import lightkurve as lk

curvy = ['TIC 166698220', 'TIC 41232189']

flat = ['TIC 229980646', 'TIC 327587572']

both_curvy_and_flat = curvy + flat

radSearch = .056 #  <-- radius in degrees = 200 arcseconds

for starName in both_curvy_and_flat:
    print('\n')
    catalogData = Catalogs.query_object(starName, radius = radSearch, catalog = "TIC")
    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']
    
    # Create a list of nearby bright stars (tess magnitude less than 14) from the rest of the data for later.
    bright = catalogData['Tmag'] < 15
    
    # Make it a list of Ra, Dec pairs of the bright ones. This is now a list of nearby bright stars.
    #nearbyStars_coordinates = list( map( lambda x,y:[x,y], catalogData[bright]['ra'], catalogData[bright]['dec'] ) )
    nearbyStars = catalogData[bright]['ID', 'Tmag', 'Jmag', 'Teff', 'logg', 'objType', 'dstArcSec']
    
    print(nearbyStars)
    
    #print('\n' + str(nearbyStars_coordinates))


#lc = lk.search_lightcurve('TIC 166698220', radius=200)

#print(lc)





















