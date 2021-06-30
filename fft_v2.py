#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 09:32:33 2021

@author: mkunz
"""
import lightkurve as lk
import scipy.fft as fft
import scipy.signal
import matplotlib.pyplot as plt
import numpy as np
#fft
lc = lk.search_lightcurve('TIC 41232189', author= 'SPOC', exptime=120)[0].download()
flat = lc.flatten().remove_nans().remove_outliers()
'''
# days to seconds
time = flat.time.value*60*60*24
flux = flat.flux.value

plt.scatter(time, flux, s=.5, c='k')
plt.xlabel('time [seconds]')
plt.ylabel('flux [e- s-]')
plt.title('Lightcurve')
plt.show()
plt.close()
'''
N=1000
a=10
time2 = np.linspace(0, a, N)
flux2 = np.sin(2*np.pi*time2/2)
x = fft.fftfreq(len(time2))[:N//2]
y = fft.fft(flux2)
plt.plot(time2, flux2)
plt.show()
plt.close()

plt.plot(x, 2/N*np.abs(y[0:N//2]))
plt.xlim(0, .01)
plt.xlabel('frequency')
plt.ylabel('power')
plt.title('FFT')
plt.show()
plt.close()
