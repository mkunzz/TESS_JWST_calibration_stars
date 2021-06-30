#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:25:10 2021

@author: mkunz
"""

import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt


tpf = lk.search_targetpixelfile("TIC 41232189", author= 'SPOC', exptime=120).download_all()

i = 0
while i < len(tpf):
    

    lc = tpf.to_lightcurve(aperture_mask='pipeline')
    lc2 = lc.remove_nans()
    lc3 = lc2.remove_outliers()
    lc4 = lc3.flatten(window_length=1501)
    
    periodogram = lc4.to_periodogram(method='lombscargle', period=np.arange(.1, 2, .001))
    best_fit_period = periodogram.period_at_max_power
    
    folded = lc4.fold(period = best_fit_period, epoch_time=.47)
    #, normalize_phase=True, epoch_time=.47
    
    print(best_fit_period)
    
    
    folded2 = folded.bin(time_bin_size=.006)
    folded2.plot()
    
    full_phase_range = folded2.phase[-1].value - folded2.phase[0].value
    tolerance = 0.05 * full_phase_range
    min_phase = folded2.time[np.argmin(folded2.flux)].value
    max_phase = folded2.time[np.argmax(folded2.flux)].value
    
    min_timestamps = folded.time_original[np.where((folded.time > min_phase - tolerance)
                                                 & (folded.time < min_phase + tolerance))].value
    max_timestamps = folded.time_original[np.where((folded.time > max_phase - tolerance)
                                                 & (folded.time < max_phase + tolerance))].value
    
    one_quarter_minima = [f for (f, t) in zip(tpf.flux.value, tpf.time.value) if t in min_timestamps]
    one_quarter_maxima = [f for (f, t) in zip(tpf.flux.value, tpf.time.value) if t in max_timestamps]
    
    avg_image = np.nanmean(tpf.flux.value, axis=0)
    diff_image = np.nanmean(one_quarter_maxima, axis=0) - np.nanmean(one_quarter_minima, axis=0)
    fig, ax = plt.subplots(1,2)
    ax[0].imshow(np.flipud(avg_image))
    ax[0].set_title('Quarter 9 average')
    ax[1].imshow(np.flipud(diff_image))
    ax[1].set_title('Quarter 9 difference image')
    fig.set_size_inches((15,6))
    
    # put colorbar
    
    
    
    
    
    
    
    
    
    
    
    
    