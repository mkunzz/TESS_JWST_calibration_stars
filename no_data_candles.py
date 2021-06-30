#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 09:07:45 2021

@author: mkunz
"""
import lightkurve as lk
#import numpy as np

no_data_std_candles = ['LDS749B', 'HD 163466', 'HD 2811', 'HD 37725', 'HD 180609 ', 'HD 55677', '1757132', '1812095 ', '1808347 ', '1802271', '1805292', '1732526 ', '1743045 ', 'HD 146233', 'HD 106252', 'HD 209458', 'HD 142331', 'P330E', 'P177D', 'C26202', 'SF1615+001A', 'SNAP-2']
'''
i=0
while i < len(no_data_std_candles):
    
    lc = lk.search_lightcurve(no_data_std_candles[i])
    print(no_data_std_candles[i])
    print(lc)
    
    i+=1
'''

lc = lk.search_lightcurve('352817378')
print('352817378')
print(lc)














                                          