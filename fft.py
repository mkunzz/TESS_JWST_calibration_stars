#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 11:29:26 2021

@author: mkunz
"""
import lightkurve as lk
import scipy.fft as fft
import scipy.signal
import matplotlib.pyplot as plt
import numpy as np

#'TIC 27533327'

lc = lk.search_lightcurve('TIC 41232189', author= 'SPOC', exptime=120)[0].download()
flat = lc.flatten().remove_nans().remove_outliers()

# days to seconds
time = flat.time.value*60*60*24
flux = flat.flux.value

plt.scatter(time, flux, s=.5, c='k')
plt.xlabel('time [seconds]')
plt.ylabel('flux [e- s-]')
plt.title('Lightcurve')
plt.show()
plt.close()

x = fft.fftfreq(len(time), 1/(time[1]-time[0]))
y = fft.fft(flux)

plt.plot(x, y)
plt.xlabel('frequency')
plt.ylabel('power')
plt.title('FFT')
plt.show()
plt.close()


#-----------------------------------------------------------------------------
#fft

def dfft(y):
    ''' based on the textbook Exercise 7.7, recursive version '''
    ''' assumes N is a power of 2 '''
    N = len(y)
    if N == 1:
        return y # (the FT of a 1-sample signal is the signal itself)
    E = dfft(y[0::2]) # even terms
    O = dfft(y[1::2]) # odd terms
    phi = np.exp(- 2j * np.pi * np.arange(N//2) / N) # twiddle factor
    return np.concatenate([E + phi*O, E - phi*O])

data = np.loadtxt("dft.txt")
tn = data[:,0]
yn = data[:,1]


plt.plot(tn,yn)
plt.title("dft.txt data")
plt.xlabel("tn")
plt.ylabel("yn")
plt.show()
plt.close()

step = tn[1]
N = len(tn)//2+1
freqs = np.linspace(0, 0.5/step, N)
powerspec = dfft(yn)
T = 2*N*step
plt.plot(freqs, (2*step**2)/T*(np.abs(powerspec)**2)[:N])
plt.show()
plt.close()

#-----------------------------------------------------------------------------

def dft(y):
    #discrete fourier transform
    N = len(y)
    c = np.zeros(N, complex)
    for k in range(N):
        for n in range(N):
            c[k] += y[n]* \
            np.exp(-2j*np.pi*k*n/N)
    return c

data = np.loadtxt("dft.txt")
tn = data[:,0]
yn = data[:,1]


plt.plot(tn,yn)
plt.title("dft.txt data")
plt.xlabel("tn")
plt.ylabel("yn")
plt.show()
plt.close()

step = tn[1]
N = len(tn)//2+1
freqs = np.linspace(0, 0.5/step, N)
powerspec = dft(yn)
T = 2*N*step

plt.plot(freqs, (2*step**2)/T*(np.abs(powerspec)**2)[:N])
'''
freqhigh = np.argmax(powerspec)
value = np.amax(powerspec)
'''
plt.title("Power spec")
plt.show()
plt.close()
pgram = scipy.signal.periodogram(yn, fs=step)

plt.plot(pgram[0], pgram[1])
plt.title("scipy periodogram")
plt.show()
plt.close()

# filter
filtered = []
for i in powerspec:
    if abs(i) < 44:
        i = 0
    else:
        pass
    filtered.append(i)
# end of filter

def inverse_dft(c):
    #discrete fourier transform
    N = len(c)
    y = np.zeros(N, complex)
    for n in range(N):
        for k in range(N):
            y[n] += c[k]* \
            np.exp(2j*np.pi*k*n/N)
    return y/N


plt.title("inverse dft back to filtered signal")
inverse = inverse_dft(filtered)
plt.plot(tn, inverse)
#plt.plot(tn, inverse_dft(powerspec))
plt.show()
plt.close()















