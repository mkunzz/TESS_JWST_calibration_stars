#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 17:07:15 2021

@author: mkunz
"""

# IMPORTS
import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from PyAstronomy.pyasl import foldAt
from scipy.optimize import curve_fit
#-----------------------------------------------------------------------------
# Search for JWST Calibration stars using
# lightkurve for Tess 2 min data using TIC numbers
lc = lk.search_lightcurve('TIC 41232189', author='SPOC', exptime=120).download_all()
#-----------------------------------------------------------------------------
# Lightcurve of single sectors
lc1_time = lc[0].time.value
lc1_flux = lc[0].flux.value
plt.scatter(lc1_time, lc1_flux, s=.6, c='k')
plt.xlabel('Time [BTJD days]')
plt.ylabel('Flux')
plt.title('Single Lightcurve from one Sector')
plt.show()
plt.close()
#%%
#-----------------------------------------------------------------------------
# Combine all sectors of data for each star
# removes nans and normalizes fluxes, adds all sector data to empty arrays
# for time and flux using 'extend'
lc2_time=[]
lc2_flux=[]
i=0
while i < len(lc):
    time_array_without_nans = lc[i].time.value[np.logical_not(np.isnan(lc[i].flux.value))]
    flux_array_without_nans = lc[i].flux.value[np.logical_not(np.isnan(lc[i].flux.value))]
    normalized_fluxes = flux_array_without_nans / np.mean(flux_array_without_nans)
    lc2_flux.extend(normalized_fluxes)
    
    lc2_time.extend(time_array_without_nans)
    
    i+=1
#%%
#-----------------------------------------------------------------------------
# Lightcurve of all sectors
plt.scatter(lc2_time, lc2_flux, s=.1, c='k')
plt.ylim(.998, 1.002)
plt.xlabel('Time [BTJD days]')
plt.ylabel('Normalized Flux')
plt.title('Combined Lightcurve from all Sectors')
plt.show()
plt.close()
#%%
#-----------------------------------------------------------------------------
# find best fit period using scipy lombscargle 
# so that lightcurve can be folded
periods = np.linspace(0.1, 1, 5000)
ang_freqs = 2 * np.pi / periods

power = signal.lombscargle(lc2_time, lc2_flux - np.mean(lc2_flux), ang_freqs, normalize=True)

N = len(lc2_time)
power *= 2 / (N * np.std(lc2_flux) ** 2)

power = signal.lombscargle(lc2_time, lc2_flux - np.mean(lc2_flux), ang_freqs)

best_fit_period = periods[np.argmax(power)]
print('Best fit period is {:.3f} days'.format(best_fit_period))
#%%
#-----------------------------------------------------------------------------
# Plot Periodograms
plt.plot(periods, power)
plt.xlabel('Period [days]')
plt.ylabel('Power')
plt.title('Lomb Scargle Periodogram\nPeriod is {:.3f} days'.format(best_fit_period))
plt.show()
plt.close()
plt.plot(ang_freqs, power)
plt.xlabel('Angular Frequency')
plt.ylabel('Power')
plt.title('Lomb Scargle Periodogram')
plt.show()
plt.close()
#%%
#-----------------------------------------------------------------------------
# Obtain the phases with respect to some
# reference point (T0)
phases = foldAt(np.array(lc2_time), best_fit_period, T0=1400.1)

# Sort with respect to phase
# First, get the order of indices ...
sortIndi = np.argsort(phases)
# ... and, second, rearrange the arrays.
phases = phases[sortIndi]
flux = np.array(lc2_flux)[sortIndi]
#%%
#-----------------------------------------------------------------------------
# Plot the folded lightcurve ie phase curve
# Overlay blue best fit curve to see variability more easily
# define the true objective function
def objective(x, a, b, c, d, e, f):
    return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f

x = phases
y = flux

# curve fit
popt, _ = curve_fit(objective, x, y)
# summarize the parameter values
a, b, c, d, e, f = popt
# plot input vs output
plt.scatter(x,y, s=.1, c='black')
# define a sequence of inputs between the smallest and largest known inputs
x_line = np.arange(min(x), max(x), .001)
# calculate the output for the range
y_line = objective(x_line, a, b, c, d, e, f)
# create a line plot for the mapping function
plt.plot(x_line, y_line, '--', color='cyan')
plt.ylim(.998, 1.002)
plt.title('Folded Lightcurve')
plt.xlabel('Phase')
plt.ylabel('Flux [e- s-1]')
plt.show()
plt.close()






















