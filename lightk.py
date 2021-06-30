#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 11:23:31 2021

@author: mkunz
"""
# import libraries
import lightkurve as lk
from lightkurve import TessTargetPixelFile
import numpy as np

#------------------------------------------------------------------------------
'''
# Target G191B2B WD
# import Tess data of G191B2B from https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
tpf = TessTargetPixelFile('/Users/mkunz/Downloads/MAST_2021-06-09T1210/TESS/tess2019331140908-s0019-0000000327587572-0164-s/tess2019331140908-s0019-0000000327587572-0164-s_tp.fits')

# aperture mask
tpf.plot(aperture_mask=tpf.pipeline_mask)

# target pixel file to lightcurve
lc = tpf.to_lightcurve()
lc.scatter()

# flatten
flat_lc = lc.flatten(window_length=1001)
flat_lc.scatter()

# periodogram to find period
periodogram = flat_lc.to_periodogram(method='bls', period=np.arange(1, 30, .001))
periodogram.plot()

best_fit_period = periodogram.period_at_max_power
print('Best fit period: {:.5f}'.format(best_fit_period))

# fold lightcurve over by best fit period
folded_lc = flat_lc.fold(period=best_fit_period)
folded_lc.scatter()

# bin to clean data
binned_lc = folded_lc.bin(binsize=10)
binned_lc.scatter()

print(lc.mission)
print(lc.targetid)

#------------------------------------------------------------------------------
# Target GD71 WD
# import Tess data of GD71 from https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
tpf = TessTargetPixelFile('/Users/mkunz/Downloads/MAST_2021-06-09T1348/TESS/tess2018349182459-s0006-0000000247923021-0126-s/tess2018349182459-s0006-0000000247923021-0126-s_tp.fits')

# aperture mask
tpf.plot(aperture_mask=tpf.pipeline_mask)

# target pixel file to lightcurve
lc = tpf.to_lightcurve()
lc.scatter()

# flatten
flat_lc = lc.flatten(window_length=1001)
flat_lc.scatter()

# periodogram to find period
periodogram = flat_lc.to_periodogram(method='bls', period=np.arange(1, 30, .001))
periodogram.plot()

best_fit_period = periodogram.period_at_max_power
print('Best fit period: {:.5f}'.format(best_fit_period))

# fold lightcurve over by best fit period
folded_lc = flat_lc.fold(period=best_fit_period)
folded_lc.scatter()

# bin to clean data
binned_lc = folded_lc.bin(binsize=10)
binned_lc.scatter()

print(lc.mission)
print(lc.targetid)

#------------------------------------------------------------------------------
# Target GD153 WD
# import Tess data of GD153 from https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
tpf = TessTargetPixelFile('/Users/mkunz/Downloads/MAST_2021-06-09T1355/TESS/tess2020078014623-s0023-0000000149505899-0177-s/tess2020078014623-s0023-0000000149505899-0177-s_tp.fits')

# aperture mask
tpf.plot(aperture_mask=tpf.pipeline_mask)

# target pixel file to lightcurve
lc = tpf.to_lightcurve()
lc.scatter()

# flatten
flat_lc = lc.flatten(window_length=1001)
flat_lc.scatter()

# periodogram to find period
periodogram = flat_lc.to_periodogram(method='bls', period=np.arange(1, 30, .001))
periodogram.plot()

best_fit_period = periodogram.period_at_max_power
print('Best fit period: {:.5f}'.format(best_fit_period))

# fold lightcurve over by best fit period
folded_lc = flat_lc.fold(period=best_fit_period)
folded_lc.scatter()

# bin to clean data
binned_lc = folded_lc.bin(binsize=10)
binned_lc.scatter()

print(lc.mission)
print(lc.targetid)

#------------------------------------------------------------------------------
# Target WD1057+719  WD
# import Tess data of WD1057+719  from https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
tpf = TessTargetPixelFile('/Users/mkunz/Downloads/MAST_2021-06-09T1401/TESS/tess2019198215352-s0014-0000000147921014-0150-s/tess2019198215352-s0014-0000000147921014-0150-s_tp.fits')

# aperture mask
tpf.plot(aperture_mask=tpf.pipeline_mask)

# target pixel file to lightcurve
lc = tpf.to_lightcurve()
lc.scatter()

# flatten
flat_lc = lc.flatten(window_length=1001)
flat_lc.scatter()

# periodogram to find period
periodogram = flat_lc.to_periodogram(method='bls', period=np.arange(1, 30, .001))
periodogram.plot()

best_fit_period = periodogram.period_at_max_power
print('Best fit period: {:.5f}'.format(best_fit_period))

# fold lightcurve over by best fit period
folded_lc = flat_lc.fold(period=best_fit_period)
folded_lc.scatter()

# bin to clean data
binned_lc = folded_lc.bin(binsize=10)
binned_lc.scatter()

print(lc.mission)
print(lc.targetid)
'''
#------------------------------------------------------------------------------
# Target WD1657+343  WD
# import Tess data of WD1657+343  from https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
tpf = TessTargetPixelFile('/Users/mkunz/Downloads/MAST_2021-06-09T1412/TESS/tess2020133194932-s0025-0000000471015233-0182-s/tess2020133194932-s0025-0000000471015233-0182-s_tp.fits')

# aperture mask
tpf.plot(aperture_mask=tpf.pipeline_mask)

# target pixel file to lightcurve
lc = tpf.to_lightcurve()
lc.scatter()

# flatten
flat_lc = lc.flatten(window_length=1001)
flat_lc.scatter()

# periodogram to find period
periodogram = flat_lc.to_periodogram(method='bls', period=np.arange(1, 30, .001))
periodogram.plot()

best_fit_period = periodogram.period_at_max_power
print('Best fit period: {:.5f}'.format(best_fit_period))

# fold lightcurve over by best fit period
folded_lc = flat_lc.fold(period=best_fit_period)
folded_lc.scatter()

# bin to clean data
binned_lc = folded_lc.bin(binsize=10)
binned_lc.scatter()

print(lc.mission)
print(lc.targetid)






















































































