#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 10:06:55 2021

@author: mkunz
"""

# find appropriate noise range

import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

lc = lk.search_lightcurve('TIC 27533327', author= 'SPOC', exptime=120)[0].download()

#-----------------------------------------------------------------------------
# Eliminate long term trends (flatten) via lightkurves's Savitzky-Golay filter
flat_lc = lc.flatten(window_length=15001)

#-----------------------------------------------------------------------------
# removes nans and normalizes fluxes

time_array_without_nans = flat_lc.time.value[np.logical_not(np.isnan(flat_lc.flux.value))]

flux_array_without_nans = flat_lc.flux.value[np.logical_not(np.isnan(flat_lc.flux.value))]
normalized_fluxes = flux_array_without_nans / np.mean(flux_array_without_nans)

#-----------------------------------------------------------------------------
# find best fit period using scipy lombscargle 
# so that lightcurve can be folded
periods = np.linspace(0.1, 5, 100000)
ang_freqs = 2 * np.pi / periods

power = signal.lombscargle(time_array_without_nans, normalized_fluxes - np.mean(normalized_fluxes), ang_freqs, normalize=False)

N = len(time_array_without_nans)
power *= 2 / (N * np.std(normalized_fluxes) ** 2)

best_fit_period = periods[np.argmax(power)]

print('Best fit period is {:.3f} days'.format(best_fit_period))

#-----------------------------------------------------------------------------
# Plot Periodograms
# Periodogram in Period Space
plt.plot(periods, power)
plt.xlabel('Period [days]')
plt.ylabel('Power')
plt.title('Lomb Scargle Periodogram\nBest Fit Period is {:.3f} days'.format(best_fit_period))
plt.show()
plt.close()

# Periodogram in Frequency Space
plt.plot(ang_freqs, power)
plt.xlabel('Angular Frequency')
plt.ylabel('Power')
plt.title('Lomb Scargle Periodogram')
plt.show()
plt.close()

#%%

print(np.std(power))

#ang_freqs[90]
#ang_freqs[2000]
































