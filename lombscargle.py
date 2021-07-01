#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 09:28:59 2021

@author: mkunz
"""

import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal


lc = lk.search_lightcurve("TIC 41232189", author= 'SPOC', exptime=120)[0].download()

lc2 = lc.flatten(window_length=1501).remove_nans().remove_outliers()
##############################################################################

# lightkurve's periodogram
pgram = lc2.to_periodogram(method='lombscargle', normalization='amplitude', maximum_period=1, oversample_factor = 10)    #maximum_frequency=120, freq_unit='microhertz',  minimum_period=.9,         

print('frequency at max amp = {:.7f}'.format(pgram.frequency_at_max_power.value))
print('max amp = {:.7f}'.format(pgram.max_power.value))
print('period at max amp = {:.7f}'.format(pgram.period_at_max_power.value) + ' days')
print('avg amp = {:.7f}'.format(np.mean(pgram.power.value)))

N = len(lc2.time.value)
sigma_rms = np.std(lc2.flux.value)
sigma_amp = np.sqrt(np.pi/N)*sigma_rms

print('sigma_amplitude = {:.7f}'.format(sigma_amp))

pgram.plot()
plt.hlines(sigma_amp, 0, 1)
plt.ylabel('Amplitude')
plt.show()
plt.close()

##############################################################################
'''
# Scipy's periodogram (scipy.signal.lombscargle)
time = lc2.time.value
flux = lc2.flux.value

# find best fit period using scipy lombscargle 
# so that lightcurve can be folded
N = len(time)
periods_1 = np.linspace(0.001, 1, N)
ang_freqs = (2*np.pi) / periods_1
mean = np.mean(flux)
pgram2 = signal.lombscargle(time, flux-mean, ang_freqs)

amplitude = np.sqrt(pgram2*(4/N))

best_fit_period = periods_1[np.argmax(pgram2)] #days
print('sector {}'.format(lc2.sector))
print('{:.5f}'.format(best_fit_period*24)+' hours')
microhertz = ang_freqs*1.8420711

plt.plot(microhertz, amplitude, c='k', lw=.5)
plt.xlim(0, 120)
plt.ylabel('Amplitude')
plt.xlabel('Frequency [uHz]')
plt.show()
plt.close()
'''
##############################################################################




