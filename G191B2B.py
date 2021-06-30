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


lk.search_lightcurve("TIC 327587572")
lc = lk.search_lightcurve("TIC 327587572")[0].download()

cleaned_lc = lc.remove_nans()
trimmed_lc = cleaned_lc.remove_outliers()

#trimmed_lc.plot()

flat_lc = trimmed_lc.flatten(window_length=1001)
#flat_lc.plot()

periodogram = flat_lc.to_periodogram(method='lombscargle', period=np.arange(.1, 10, .001))
#periodogram.plot()

best_fit_period = periodogram.period_at_max_power
print('Best fit period: {:.5f}'.format(best_fit_period))

# fold lightcurve over by best fit period
folded_lc = flat_lc.fold(period=best_fit_period)

#folded_lc.scatter()
#plt.savefig('/Users/mkunz/figures/' + 'test_star.png')
'''
# bin to clean data
binned_lc = folded_lc.bin(binsize=10)
binned_lc.scatter()
'''
print(lc.mission)
print(lc.targetid)

#print(np.std(folded_lc.flux))
x = folded_lc.time.value
y = folded_lc.flux
 
# define the true objective function
def objective(x, a, b, c, d, e, f):
	return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f
 
# curve fit
popt, _ = curve_fit(objective, x, y)
# summarize the parameter values
a, b, c, d, e, f = popt
# plot input vs output
plt.scatter(x,y, s=.2, c='black')
# define a sequence of inputs between the smallest and largest known inputs
x_line = np.arange(min(x), max(x), .001)
# calculate the output for the range
y_line = objective(x_line, a, b, c, d, e, f)
# create a line plot for the mapping function
plt.plot(x_line, y_line, '--', color='red')










