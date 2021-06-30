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
pgram = lc2.to_periodogram(method='lombscargle', normalization='amplitude', oversample_factor = 10, maximum_frequency=120, freq_unit='microhertz')
pgram.plot()
plt.ylabel('Amplitude')
plt.show()
plt.close()
##############################################################################
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

#best_fit_period = periods_1[np.argmax(pgram2)]

microhertz = ang_freqs*1.8420711

plt.plot(microhertz, amplitude, c='k', lw=.5)
plt.xlim(0, 120)
plt.ylabel('Amplitude')
plt.xlabel('Frequency [uHz]')
plt.show()
plt.close()