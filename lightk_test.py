# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#imports
import lightkurve as lk
import pandas
import numpy as np
import astropy
import astroquery
import matplotlib.pyplot as plt

lk.search_lightcurve("Kepler-10", mission="Kepler", quarter = 6)
lc = lk.search_lightcurve("Kepler-10", mission="Kepler", quarter = 6)[3].download()
lc.plot()

##########

from lightkurve import search_targetpixelfile
pixelfile = search_targetpixelfile('KIC 8462852', quarter=16).download(quality_bitmask='hardest')
pixelfile.plot(frame=1)
lc = pixelfile.to_lightcurve(aperture_mask='all')

# object lc; lc.time , lc.flux

lc.plot()

pixelFile = search_targetpixelfile('KIC 6922244', quarter=4).download()
lc = pixelFile.to_lightcurve(aperture_mask=pixelFile.pipeline_mask)
lc.plot()

flat_lc = lc.flatten(window_length=401)
flat_lc.plot()

folded_lc = flat_lc.fold(period=3.5225)
folded_lc.plot()

# period somewhere between .3 to 10 days
periodogram = lc.to_periodogram(method='bls', period=np.arange(.3, 10, .001))
periodogram.plot()

best_fit_period = periodogram.period_at_max_power
print('Best fit period: {:.5f}'.format(best_fit_period))

binned_lc = folded_lc.bin(binsize=10)
binned_lc.plot()

lc.remove_nans().flatten(window_length=401).fold(period=3.5225).bin(binsize=10).plot()

###########

from lightkurve import search_targetpixelfile
import lightkurve as lk
search_result = lk.search_targetpixelfile('Pi Mensae', mission='TESS', sector=1)
tpf = search_result.download(quality_bitmask='default')

# tpf.mission , tpf.targetid

aperture_mask = tpf.create_threshold_mask(threshold=10)
lc = tpf.to_lightcurve(aperture_mask=aperture_mask)
lc.scatter()

flat_lc = lc.flatten(window_length=1001)
flat_lc.errorbar()

# remove telescope moving data

mask = (flat_lc.time < 1346) | (flat_lc.time > 1350)
masked_lc = flat_lc[mask]
masked_lc.errorbar()

masked_lc.scatter(s=0.1)

# remove outliers
clipped_lc = masked_lc.remove_outliers(sigma=6)
clipped_lc.errorbar()

clipped_lc.scatter(s=0.1)

########

from lightkurve import TessTargetPixelFile
import lightkurve as lk
#from lightkurve import TessLightCurveFile
tpf = TessTargetPixelFile('....... lc.fits')

















