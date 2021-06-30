#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 10:02:15 2021

@author: mkunz
"""

import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

print('Length of candles with TIC = ' + str(len(candles_with_TIC)))
print('Length of candles with 2 minute data = ' + str(len(candles_with_2min)))

# define the true objective function
def objective(x, a, b, c, d, e, f):
    return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f

iii = 0
while iii < len(candles_with_2min):
    
    lc = lk.search_lightcurve(candles_with_2min[iii], author='SPOC', exptime=120).download_all()

    ii=0
    while ii < len(lc):
        cleaned_lc = lc[ii].remove_nans()
        flat_lc = cleaned_lc.flatten(window_length=1501)
        trimmed_lc = flat_lc.remove_outliers()
        
        trimmed_lc.scatter()
        plt.title('Lightcurve of ' + str(candles_with_2min[iii]) + '\n' + 'Sector ' + str(lc[ii].sector))
        plt.savefig('/Users/mkunz/figures/lc' + str(candles_with_2min[iii]) + 'sector' + str(lc[ii].sector) + '.png')
        plt.show()
        plt.close()
        
        periodogram = trimmed_lc.to_periodogram(method='lombscargle', frequency=np.arange(0.1, 10, .001))
        best_fit_period = periodogram.period_at_max_power
        periodogram.plot()
        plt.title('Periodogram of ' + str(candles_with_2min[iii])  + '\n' + 'Sector ' + str(lc[ii].sector) + '; Period = {:.3f}'.format(best_fit_period))
        plt.savefig('/Users/mkunz/figures/periodogram' + str(candles_with_2min[iii]) + 'sector' + str(lc[ii].sector) + '.png')
        plt.show()
        plt.close()
        
        folded_lc = trimmed_lc.fold(period=best_fit_period, normalize_phase=True)
        
        #print(np.std(folded_lc.flux))
        
        x = folded_lc.time.value
        y = folded_lc.flux
         
        # curve fit
        popt, _ = curve_fit(objective, x, y)
        # summarize the parameter values
        a, b, c, d, e, f = popt
        # plot input vs output
        plt.scatter(x,y, s=.4, c='black')
        # define a sequence of inputs between the smallest and largest known inputs
        x_line = np.arange(min(x), max(x), .001)
        # calculate the output for the range
        y_line = objective(x_line, a, b, c, d, e, f)
        # create a line plot for the mapping function
        plt.plot(x_line, y_line, '--', color='cyan')
        plt.xlabel('Phase [JD]')
        plt.ylabel('Flux [e- s-1]')
        plt.title('Star : ' + str(candles_with_2min[iii]) + '\n' +  'Sector: ' + str(lc[ii].sector) + ' ; Period: {:.3f}'.format(best_fit_period))
        plt.savefig('/Users/mkunz/figures/' + str(candles_with_2min[iii]) + 'sector' + str(lc[ii].sector) + '.png')
        plt.show()
        plt.close()
        
        ii+=1

    iii+=1

























