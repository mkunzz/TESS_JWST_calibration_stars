#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 09:28:59 2021

@author: mkunz
"""

import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt


tpf = lk.search_targetpixelfile("TIC 41232189", author= 'SPOC', exptime=120)[1].download()

def sin:
    # fit sin wave
    # epoch can vary
    # fix period, x is time, y is flux
    return A * sine( 2 * pi * period / time + epoch)

'''
#tpf.plot()
index = 0
powers= []
i = 0
while i < 13:
    j = 0
    while j < 13:
        # Choosing specific pixels to use by defining a new aperture mask manually
        aper_mask = np.zeros(tpf.shape[1:], dtype=bool)
        aper_mask[5+i:8+i, 5+j:8+j] = True
        tpf.plot(aperture_mask=aper_mask, mask_color='red')
        plt.show()
        plt.close()
        
        
        lc = tpf.to_lightcurve(aperture_mask=aper_mask)
        lc_new = lc.remove_nans()
        flat_lc = lc_new.flatten(window_length=1501)
        trimmed_lc = flat_lc.remove_outliers()
        
        #trimmed_lc.scatter()
        
        
        periodogram = trimmed_lc.to_periodogram(method='lombscargle', period=np.arange(.1, 2, .001))
        best_fit_period = periodogram.period_at_max_power
        #periodogram.plot()
        #plt.show()
        #plt.close()
        print('i = ' +str(i))
        print('j = ' + str(j))
        print('index is ' + str(index))
        print(periodogram.max_power)
        powers.append(periodogram.max_power.value)
        
        index += 1
        
        j+=1
    i+=1

print('the max power is ' + str(np.max(powers)))
print('the index of max power is ' + str(np.argmax(powers)))
'''


# Choosing specific pixels to use by defining a new aperture mask manually
aper_mask = np.zeros(tpf.shape[1:], dtype=bool)
aper_mask[12:15, 17:20] = True
tpf.plot(aperture_mask=aper_mask, mask_color='red')



lc = tpf.to_lightcurve(aperture_mask=aper_mask)
lc_new = lc.remove_nans()
flat_lc = lc_new.flatten(window_length=1501)
trimmed_lc = flat_lc.remove_outliers()

trimmed_lc.scatter()


periodogram = trimmed_lc.to_periodogram(method='lombscargle', period=np.arange(.1, 2, .001))
best_fit_period = periodogram.period_at_max_power
periodogram.plot()
print(best_fit_period)

print(periodogram.max_power)

folded = trimmed_lc.fold(period=best_fit_period, normalize_phase=True)
folded.scatter()

#periodogram of each pixel
#tpf.plot_pixels(aperture_mask='pipeline', periodogram=True)



