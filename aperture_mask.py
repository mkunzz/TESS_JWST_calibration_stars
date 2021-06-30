#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 14:37:07 2021

@author: mkunz
"""

import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


tpf = lk.search_targetpixelfile("TIC 41232189", author= 'SPOC', exptime=120)[1].download()

tpf.plot()
'''
# Defining the aperture mask using the threshold method
target_mask = tpf.create_threshold_mask(threshold=800, reference_pixel='center')
n_target_pixels = target_mask.sum()
print(n_target_pixels)
'''
tpf.plot(aperture_mask='pipeline', mask_color='r')

lc = tpf.to_lightcurve(aperture_mask='pipeline')
lc = lc.remove_nans()
lc = lc.remove_outliers()
target_lc = lc.flatten(window_length=1501)
target_lc.scatter()

periodogram = target_lc.to_periodogram(method='lombscargle', period=np.arange(.1, 2, .001))
best_fit_period = periodogram.period_at_max_power
periodogram.plot()
print(best_fit_period)
print(periodogram.max_power)

folded = target_lc.fold(period = best_fit_period, normalize_phase=True)
folded.scatter()

'''
# Choosing specific pixels to use by defining a new aperture mask manually
aper_new = np.zeros(tpf.shape[1:], dtype=bool)
aper_new[6:10, 1:4] = True
tpf.plot(aperture_mask=aper_new, mask_color='red')



#target_lc = tpf.to_lightcurve(aperture_mask=target_mask)
target_lc = tpf.to_lightcurve(aperture_mask='threshold')

target_lc = target_lc.flatten(window_length=1501)
target_lc = target_lc.remove_outliers()
target_lc = target_lc.remove_nans()
target_lc.scatter(label='labels...')

periodogram = target_lc.to_periodogram(method='lombscargle', period=np.arange(.1, 10, .001))
best_fit_period = periodogram.period_at_max_power
periodogram.plot()
print(best_fit_period)

folded_lc = target_lc.fold(period=best_fit_period, normalize_phase=True)

folded_lc.scatter()

##############################################################################

lc = lk.search_lightcurve("TIC 41232189", author= 'SPOC', exptime=120)[2].download()

lc = lc.flatten(window_length=1001)
lc = lc.remove_outliers()
lc = lc.remove_nans()
lc.scatter()

periodogram = lc.to_periodogram(method='lombscargle', period=np.arange(.1, 10, .001))
best_fit_period = periodogram.period_at_max_power
periodogram.plot()
print(best_fit_period)
folded = lc.fold(period=best_fit_period, normalize_phase=True)
folded.scatter()
'''




