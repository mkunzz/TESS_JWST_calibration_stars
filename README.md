# TESS_JWST_calibration_stars
Using TESS to analyze the variability of the JWST candidate calibration stars.  

# All of the Functions to analyze JWST Calibration Stars.

    They are all inside of michael.py

Example: 

    [In]: import michael

    [In]: michael.get_sectors(41232189)
    
    [Out]: [3, 6, 9, 10, 13, 30, 33, 36, 37]

This returned all of the TESS sectors available for the JWST standard star TIC 41232189.

# These functions will:

Function | Description
------------ | -------------
get_stars(data_type) | lists all stars of desired type; '2min', '30min', 'variables', 'both', 'all', 'none'
get_sectors(tic_number) | get all available sectors for 1 star
lc_1sector(tic_number, sector_number) | return light curve object for 1star, 1sector
plot_lightcurve(tic_number, sector_number) | plot lightcurve for 1star, 1sector
plot_all_lightcurves() | return lightcurve plots of all stars with TESS data
upper_limit_percent_variation(tic_number, sector_number) | gives upper limit on variability in percent (for variable stars); max - min flux
flux_variation(tic_number, sector_number) | flux variation in percent via 99 percentile - 1 percentile
pgram(tic_number, sector_number) | returns max amp [%] and best period [hrs]
plot_periodogram(tic_number, sector_number, max_freq) | plots the periodogram in frequency in uHz and period in hours. asks for max freq (x limit on plot)
stitched_pgram(tic_number, min_freq, max_freq) | plot stitched pgram of all sectors
difference_imaging(tic_number, sector_number, best_period, epoch_time, binsize, tolerance, avg_vmax, diff_vmax) | perform difference imaging
NearbyBrightStars(tic_number) | list nearby bright stars < 15 Tmag
saving_eleanor_data_to_fits_file_locally(tic_number) | save eleanor 30min data to fits files locally
see_background_stars_in_tpf(tic_number, sector_number) | make an interactive plot of tpf with background stars plotted (with their stats) as well
two_or_thirty_or_none(tic_number) | tells you if the star has 2min data or 30min data or none
info_on_single_star(tic_number) | give back information about star with input tic number; TIC, ra, dec, Tmag, Vmag, Kmag, Teff, log g
info_on_all_stars(name_of_csv) | info on all stars in csv (excel format) -> TIC, ra, dec, Tmag, Vmag, Kmag, Teff, log g
fold_lc(tic_number, sector_number) | plots folded lightcurve
crowd_sap_1star_1sector(tic_number, sector_number) | returns the crowd sap value for 1 sector of 1 star
crowd_sap_all_2min(array) | crowd sap (for all 2min data), writes to a csv (excel format)
bootstrap(tic_number, sector_number) | bootstrap method to give 5*avg noise level for 1 sector of 1 star
stitch(tic_number) | stitches together and plots all sectors for 1 star
one_noise_level(N, fluxes)| equation to return 1*avg noise level; as inputs it needs N number of data points in time series (ie len(lc.time.value)) and the flux array (ie lc.flux.value)
bootstrap_all_standards(name_of_csv) | bootstrap method to give 5*avg noise level for all JWST standard stars
std_or_diff_for_variables(tic_number) | gives back std of max amplitudes [%] for all sectors, std of period [hours] for all sectors (if 3 or more sectors std, if 2 sectors calculates difference not std)

