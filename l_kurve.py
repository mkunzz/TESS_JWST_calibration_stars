#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 14:58:39 2021

@author: mkunz
"""

import lightkurve as lk
import numpy as np

WD_standard_candles = ['G191B2B', 'GD71', 'GD153 ', 'LDS749B', 'WD1057+719', 'WD1657+343']
A_star_standard_candles = ['HD 166205', 'HD 14943', 'HD 42525', 'HD 128998', 'HD 158485', 'HD 163466', 'HD 101452', 'HD 2811', 'HD 37725', 'HD 116405', 'HD 180609', 'HD 55677', 'BD+60 1753', '1757132', '1812095', '1808347', '1802271', '1805292', '1732526', '1743045']
G_star_standard_candles = ['HD 146233', 'HD 186427', 'HD 159222', 'HD 205905', 'HD 106252', 'HD 37962', 'HD 209458', 'HD 38949', 'HD 142331', 'HD 167060', 'HD 115169', 'P330E', 'P177D', 'C26202', 'SF1615+001A', 'SNAP-2']
all_std_candles = WD_standard_candles+A_star_standard_candles+G_star_standard_candles


'''
i = 0
while i < 42:
    print('-------------------------------------------------------------------------------')
    print('For the standard candle star: ' + all_std_candles[i])
    print(lk.search_lightcurve(all_std_candles[i], exptime=120))
    i+=1
'''

std_candles_with_2min_data = ['G191B2B', 'GD71', 'GD153 ', 'WD1057+719', 'WD1657+343', 'HD 166205', 'HD 14943', 'HD 42525', 'HD 128998', 'HD 158485', 'HD 101452', 'HD 116405', 'BD+60 1753', 'HD 186427', 'HD 159222', 'HD 205905', 'HD 37962', 'HD 38949', 'HD 167060', 'HD 115169']

i = 0
while i < len(std_candles_with_2min_data):
    print('-------------------------------------------------------------------------------')
    print('For the standard candle star: ' + std_candles_with_2min_data[i])
    lc = lk.search_lightcurve(std_candles_with_2min_data[i])[0].download()
    #lc.plot()
    
    flat_lc = lc.flatten(window_length=1001)
    #flat_lc.plot()
    
    periodogram = flat_lc.to_periodogram(method='lombscargle', period=np.arange(1, 10, .001))
    #periodogram.plot()
    
    best_fit_period = periodogram.period_at_max_power
    print('Best fit period: {:.5f}'.format(best_fit_period))
    
    # fold lightcurve over by best fit period
    folded_lc = flat_lc.fold(period=best_fit_period)
    #folded_lc.scatter()
    
    # bin to clean data
    binned_lc = folded_lc.bin(binsize=10)
    binned_lc.scatter()
    
    print(lc.mission)
    print(lc.targetid)

    i+=1














