#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 16:37:57 2021

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
# define the true objective function
def objective(x, a, b, c, d, e, f):
    return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f
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

#print('Length of candles with TIC = ' + str(len(candles_with_TIC)))
#print('Length of candles with 2 minute data = ' + str(len(candles_with_2min)))
#print('-----------------------------------------------------------------------')
#-----------------------------------------------------------------------------
ab = 0
while ab < len(candles_with_2min):
    
    #print('Loop number {}'.format(ab+1))
    #print('out of {}'.format(len(candles_with_2min)))
    
    lc = lk.search_lightcurve(candles_with_2min[ab], author='SPOC', exptime=120).download_all()
    '''
    #-----------------------------------------------------------------------------
    # Eliminate long term trends (flatten) via Savitzky-Golay filter
    flat_lc = lc[0].flatten(window_length=1501)
    #-----------------------------------------------------------------------------
    
    # Lightcurve of single sector
    lc1_time = flat_lc.time.value
    lc1_flux = flat_lc.flux.value
    plt.scatter(lc1_time, lc1_flux, s=.9, c='k')
    plt.xlabel('Time [BTJD days]')
    plt.ylabel('Flux [e- s-]')
    plt.title('Lightcurve of Single Sector\nof Star {}'.format(candles_with_2min[ab]))
    plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[ab]) + 'single_lc' + '.png')
    plt.show()
    plt.close()
    '''
    #-----------------------------------------------------------------------------
    # Combine all sectors of data for each star
    # removes nans and normalizes fluxes, adds all sector data to empty arrays
    # for time and flux using 'extend'
    lc2_time=[]
    lc2_flux=[]
    n=0
    while n < len(lc):
        
        flat_lc = lc[n].flatten(window_length=1501).remove_outliers()
        
        time_array_without_nans = flat_lc.time.value[np.logical_not(np.isnan(flat_lc.flux.value))]
        
        flux_array_without_nans = flat_lc.flux.value[np.logical_not(np.isnan(flat_lc.flux.value))]
        normalized_fluxes = flux_array_without_nans / np.mean(flux_array_without_nans)
        
        lc2_flux.extend(normalized_fluxes)
        lc2_time.extend(time_array_without_nans)
        
        n+=1
    '''
    #-----------------------------------------------------------------------------
    # Lightcurve of all sectors
    plt.scatter(lc2_time, lc2_flux, s=.9, c='k', label=candles_with_2min[ab])
    plt.ylim(np.mean(lc2_flux)-3*np.std(lc2_flux), np.mean(lc2_flux)+3*np.std(lc2_flux))
    plt.xlabel('Time [BTJD days]')
    plt.ylabel('Normalized Flux [e- s-]')
    plt.title('Combined Lightcurve from all Sectors\nNumber of Sectors = {}'.format(len(lc)))
    plt.legend()
    plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[ab]) + 'combined_lc' + '.png')
    plt.show()
    plt.close()
    '''
    #-----------------------------------------------------------------------------
    print('For star {}, the stdev of lc is '.format(candles_with_2min[ab]) + str(np.std(lc2_flux)))
    #-----------------------------------------------------------------------------
    # find best fit period using scipy lombscargle 
    # so that lightcurve can be folded
    periods_1 = np.linspace(0.001, 1, 100000)
    ang_freqs = 2 * np.pi / periods_1
    
    power = signal.lombscargle(lc2_time, lc2_flux - np.mean(lc2_flux), ang_freqs, normalize=False)
    
    N = len(lc2_time)
    power *= 2 / (N * np.std(lc2_flux) ** 2)
    
    best_fit_period = periods_1[np.argmax(power)]
    #print('{}'.format(candles_with_2min[ab]))
    print('Best fit period is {:.3f} days [range .001-1 days]'.format(best_fit_period))
    print('average power: ' + str(np.mean(power)))
    print('max power: ' + str(np.max(power)))
    '''
    #-----------------------------------------------------------------------------
    # Plot Periodograms
    # Periodogram in Period Space
    plt.plot(periods_1, power, label=candles_with_2min[ab])
    plt.xlabel('Period [days]')
    plt.ylabel('Power')
    plt.title('Lomb Scargle Periodogram\nBest Fit Period is {:.3f} days [range .001-1 days]'.format(best_fit_period))
    plt.legend()
    plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[ab]) + 'period_pgram[range .001-1 days]' + '.png')
    plt.show()
    plt.close()
    # Periodogram in Frequency Space
    plt.plot(ang_freqs, power)
    plt.xlabel('Angular Frequency')
    plt.ylabel('Power')
    plt.title('Lomb Scargle Periodogram\nof Star {} \n[range .001-1 days]'.format(candles_with_2min[ab]))
    plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[ab]) + 'freq_pgram[range .001-1 days]' + '.png')
    plt.show()
    plt.close()
    #-----------------------------------------------------------------------------
    # Obtain the phases with respect to some
    # reference point (T0)
    phases = foldAt(np.array(lc2_time), best_fit_period, T0=np.mean(lc2_time))
    
    # Sort with respect to phase
    # First, get the order of indices ...
    sortIndi = np.argsort(phases)
    # ... and, second, rearrange the arrays.
    phases = phases[sortIndi]
    flux = np.array(lc2_flux)[sortIndi]
    
    #-----------------------------------------------------------------------------
    # Plot the folded lightcurve ie phase curve
    # Overlay blue best fit curve to see variability more easily
   
    x = phases
    y = flux
    
    # curve fit
    popt, _ = curve_fit(objective, x, y)
    # summarize the parameter values
    a, b, c, d, e, f = popt
    # plot input vs output
    plt.scatter(x,y, s=.9, c='black')
    # define a sequence of inputs between the smallest and largest known inputs
    x_line = np.arange(min(x), max(x), .001)
    # calculate the output for the range
    y_line = objective(x_line, a, b, c, d, e, f)
    
    # create a line plot for the mapping function
    plt.plot(x_line, y_line, '--', color='red')
    plt.ylim(np.mean(lc2_flux)-3*np.std(lc2_flux), np.mean(lc2_flux)+3*np.std(lc2_flux))
    plt.title('Folded Lightcurve\nof Star {} \n[range .001-1 days]'.format(candles_with_2min[ab]))
    plt.xlabel('Phase')
    plt.ylabel('Flux [e- s-1]')
    plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[ab]) + 'folded_lc[range .001-1 days]' + '.png')
    plt.show()
    plt.close()
    '''
    print('------------------------------------------------------------------')
    periods_2 = np.linspace(0.9, 5, 100000)
    ang_freqs = 2 * np.pi / periods_2
    
    power = signal.lombscargle(lc2_time, lc2_flux - np.mean(lc2_flux), ang_freqs, normalize=False)
    
    N = len(lc2_time)
    power *= 2 / (N * np.std(lc2_flux) ** 2)
    
    best_fit_period = periods_2[np.argmax(power)]
    print('{}'.format(candles_with_2min[ab]))
    print('Best fit period is {:.3f} days [range .9-5 days]'.format(best_fit_period))
    print('average power: ' + str(np.mean(power)))
    print('max power: ' + str(np.max(power)))
    #-----------------------------------------------------------------------------
    '''
    # Plot Periodograms
    # Periodogram in Period Space
    plt.plot(periods_2, power, label=candles_with_2min[ab])
    plt.xlabel('Period [days]')
    plt.ylabel('Power')
    plt.title('Lomb Scargle Periodogram\nBest Fit Period is {:.3f} days [range .9-5 days]'.format(best_fit_period))
    plt.legend()
    plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[ab]) + 'period_pgram[range .9-5 days]' + '.png')
    plt.show()
    plt.close()
    # Periodogram in Frequency Space
    plt.plot(ang_freqs, power)
    plt.xlabel('Angular Frequency')
    plt.ylabel('Power')
    plt.title('Lomb Scargle Periodogram\nof Star {} \n[range .9-5 days]'.format(candles_with_2min[ab]))
    plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[ab]) + 'freq_pgram[range .9-5 days]' + '.png')
    plt.show()
    plt.close()
    #-----------------------------------------------------------------------------
    # Obtain the phases with respect to some
    # reference point (T0)
    phases = foldAt(np.array(lc2_time), best_fit_period, T0=np.mean(lc2_time))
    
    # Sort with respect to phase
    # First, get the order of indices ...
    sortIndi = np.argsort(phases)
    # ... and, second, rearrange the arrays.
    phases = phases[sortIndi]
    flux = np.array(lc2_flux)[sortIndi]
    
    #-----------------------------------------------------------------------------
    # Plot the folded lightcurve ie phase curve
    # Overlay blue best fit curve to see variability more easily
    
    x = phases
    y = flux
    
    # curve fit
    popt, _ = curve_fit(objective, x, y)
    # summarize the parameter values
    a, b, c, d, e, f = popt
    # plot input vs output
    plt.scatter(x,y, s=.9, c='black')
    # define a sequence of inputs between the smallest and largest known inputs
    x_line = np.arange(min(x), max(x), .001)
    # calculate the output for the range
    y_line = objective(x_line, a, b, c, d, e, f)
    
    # create a line plot for the mapping function
    plt.plot(x_line, y_line, '--', color='red')
    plt.ylim(np.mean(lc2_flux)-3*np.std(lc2_flux), np.mean(lc2_flux)+3*np.std(lc2_flux))
    plt.title('Folded Lightcurve\nof Star {} \n[range .9-5 days]'.format(candles_with_2min[ab]))
    plt.xlabel('Phase')
    plt.ylabel('Flux [e- s-1]')
    plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[ab]) + 'folded_lc[range .9-5 days]' + '.png')
    plt.show()
    plt.close()
    '''
    print('------------------------------------------------------------------')
    
    ab+=1




       

