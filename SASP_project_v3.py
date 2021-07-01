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
#import matplotlib.pyplot as plt
#import scipy.signal as signal
#-----------------------------------------------------------------------------

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
variables = []
with open('2mindata.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',', dialect='excel')
    spamwriter.writerow(['star', 'sector', 'std of lc', 'best period [days] [range 0-1 days]', 'max amplitude [e- s-] [range 0-1 days]', 'average amplitude [range 0-1 days]', 'best period [days] [range .9-5 days]', 'max amplitude [e- s-] [range .9-5 days]', 'average amplitude [range .9-5 days]', 'mean noise level [sigma_amp]'])
    ab = 0
    while ab < len(candles_with_2min):
        
        lc = lk.search_lightcurve(candles_with_2min[ab], author='SPOC', exptime=120).download_all()
    
        n=0
        while n < len(lc):
            
            flat_lc = lc[n].flatten(window_length=1501).remove_outliers().remove_nans()
            # lightkurve's periodogram
            # periodogram 1 from 0-1 days
            # periodogram 2 from .9-5 days
            pgram = flat_lc.to_periodogram(method='lombscargle', normalization='amplitude', maximum_period=1, oversample_factor = 10)  #freq_unit='microhertz'
            pgram2 = flat_lc.to_periodogram(method='lombscargle', normalization='amplitude', minimum_period=.9, maximum_period=5, oversample_factor = 10)  #freq_unit='microhertz'
            
            N = len(flat_lc.time.value)
            sigma_rms = np.std(flat_lc.flux.value)
            # mean noise level in amplitude spectrum
            sigma_amp = np.sqrt(np.pi/N)*sigma_rms
            '''
            # plotting periodograms for 0-1 days and 1-5 days with avg noise level in blue
            pgram.plot()
            plt.hlines(sigma_amp, 0, 1)
            plt.ylabel('Amplitude')
            plt.show()
            plt.close()
            
            pgram2.plot()
            plt.hlines(sigma_amp, .9, 5)
            plt.ylabel('Amplitude')
            plt.show()
            plt.close()
            '''
            a = candles_with_2min[ab]
            b = flat_lc.sector
            c = np.std(flat_lc.flux.value)
            d = pgram.period_at_max_power.value
            e = pgram.max_power.value
            f = np.mean(pgram.power.value)
            g = pgram2.period_at_max_power.value
            h = pgram2.max_power.value
            i = np.mean(pgram2.power.value)
            j = sigma_amp
            
            spamwriter.writerow([a] + [b] + [c] + [d] + [e] + [f] + [g] + [h] + [i] + [j])
            
            if e > 3*j:
                variables.append(a)
                #print('star ' + a + ', sector ' + str(b) + ' range 0-1 days, detected variable: max > 3 sigma')
            else:
                pass
            if h > 3*j:
                variables.append(a)
                #print('star ' + a + ', sector ' + str(b) + ' range 1-5 days, detected variable: max > 3 sigma')
            else:
                pass
            
            n+=1
        
        ab+=1

variables2 = []
for star in variables:
    if star not in variables2:
        variables2.append(star)
print(len(variables2))































