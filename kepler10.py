#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 12:26:40 2021

@author: mkunz
"""
import lightkurve as lk
import numpy as np

lk.search_lightcurve("Kepler-10", mission="Kepler", quarter = 6)
lc = lk.search_lightcurve("Kepler-10", mission="Kepler", quarter = 6)[3].download()
lc.plot()

# flatten
flat_lc = lc.flatten(window_length=1001)
flat_lc.scatter()

# periodogram to find period
periodogram = flat_lc.to_periodogram(method='bls', period=np.arange(1, 20, .001))
periodogram.plot()

best_fit_period = periodogram.period_at_max_power
print('Best fit period: {:.5f}'.format(best_fit_period))

# fold lightcurve over by best fit period
folded_lc = flat_lc.fold(period=best_fit_period)
folded_lc.scatter()
'''
# bin to clean data
binned_lc = folded_lc.bin(binsize=10)
binned_lc.scatter()
'''

print(lc.mission)
print(lc.targetid)