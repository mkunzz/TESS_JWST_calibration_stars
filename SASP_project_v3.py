#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:40:51 2021

@author: mkunz
"""
# IMPORTS
import csv
import lightkurve as lk
import numpy as np

#-----------------------------------------------------------------------------
# define the true objective function
def objective(x, a, b, c, d, e, f):
    return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f
#-----------------------------------------------------------------------------
with open('2mindata.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    i = 0
    while i < 10:
        spamwriter.writerow([i] * 5 + ['Baked Beans'])
        spamwriter.writerow(['Spam', 'Lovely Spam', 'Wonderful Spam'])
        
        i+=1
#-----------------------------------------------------------------------------
# Search for JWST Calibration stars using
# lightkurve for Tess 2 min data using TIC numbers

candles_with_TIC = ['TIC 327587572', 'TIC 247923021', 'TIC 149505899', 'TIC 352817378', 'TIC 147921014', 'TIC 471015233', 'TIC 383553764', 'TIC 166698220', 'TIC 41232189', 'TIC 298165335', 'TIC 198456033', 'TIC 219820925', 'TIC 181240911', 'TIC 80313923', 'TIC 75586606', 'TIC 165370459', 'TIC 229945862', 'TIC 440765193', 'TIC 219752116', 'TIC 135656809', 'TIC 27533327', 'TIC 39464221', 'TIC 441120034', 'TIC 397558558', 'TIC 140282069', 'TIC 420814525', 'TIC 32869782', 'TIC 54837036', 'TIC 365653206', 'TIC 229980646', 'TIC 8591766', 'TIC 417544924', 'TIC 144599609', 'TIC 315627636', 'TIC 207440438', 'TIC 219094190', 'TIC 233095291', 'TIC 219114641', 'TIC 233067231', 'TIC 233075513', 'TIC 219897252', 'TIC 233205654']
candles_with_2min = []
i = 0
while i < len(candles_with_TIC):
    lc = lk.search_lightcurve(candles_with_TIC[i], author='SPOC', exptime=120)
    if len(lc) == 0:
        pass
    else:
        candles_with_2min.append(candles_with_TIC[i])
    i+=1
#-----------------------------------------------------------------------------
ab = 0
while ab < len(candles_with_2min):
    
    #print('Loop number {}'.format(ab+1))
    #print('out of {}'.format(len(candles_with_2min)))
    
    lc = lk.search_lightcurve(candles_with_2min[ab], author='SPOC', exptime=120).download_all()
  
    #-----------------------------------------------------------------------------
    n=0
    while n < len(lc):
        
        flat_lc = lc[n].flatten(window_length=1501).remove_outliers().remove_nans()
        flat_lc.periodogram()
        n+=1
    
    ab+=1

