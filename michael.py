##################################################################################
# Imports
import eleanor
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from astroquery.mast import Catalogs
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
import csv
import random
####################################################################################
# Arrays of all of the JWST standards
two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
              41232189,298165335,198456033,181240911,165370459,219752116,27533327,
              39464221,441120034,140282069,32869782,365653206,229980646]
thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                 219094190,233095291,219114641,233067231,233075513,219897252,233205654]
none = [352817378,135656809,397558558,420814525,54837036,315627636,80313923,75586606]
variable_list = [383553764,166698220,41232189,198456033,441120034,32869782,219820925]
full_list_of_JWST_standards = two_minute+thirty_minute+none
##################################################################################
# Functions
def get_stars(data_type):
    array = []
    if data_type == '2min':
        array.append([327587572,247923021,149505899,147921014,471015233,383553764,166698220,
              41232189,298165335,198456033,181240911,165370459,219752116,27533327,
              39464221,441120034,140282069,32869782,365653206,229980646])
    elif data_type == '30min':
        array.append([219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                 219094190,233095291,219114641,233067231,233075513,219897252,233205654])
    elif data_type == 'variables':
        array.append([383553764,166698220,41232189,198456033,441120034,32869782,219820925])
    elif data_type == 'both':
        array.append([327587572,247923021,149505899,147921014,471015233,383553764,166698220,
              41232189,298165335,198456033,181240911,165370459,219752116,27533327,
              39464221,441120034,140282069,32869782,365653206,229980646,219820925,229945862,
              440765193,8591766,417544924,144599609,207440438,219094190,233095291,219114641,
              233067231,233075513,219897252,233205654])
    elif data_type == 'none':
        array.append([352817378,135656809,397558558,420814525,54837036,315627636,80313923,75586606])
    elif data_type == 'all':
        array.append([327587572,247923021,149505899,147921014,471015233,383553764,166698220,
              41232189,298165335,198456033,181240911,165370459,219752116,27533327,
              39464221,441120034,140282069,32869782,365653206,229980646,219820925,229945862,
              440765193,8591766,417544924,144599609,207440438,219094190,233095291,219114641,
              233067231,233075513,219897252,233205654,352817378,135656809,397558558,
              420814525,54837036,315627636,80313923,75586606])
    else:
        pass
    return array[0]
        
def get_sectors(tic_number):
    # Get available sectors for JWST calibration stars.
    # Works for 2min and 30min stars
    # Input is TIC number.
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    
    sectors_all = []
    
    if tic_number in two_minute:
        lc_all = lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC', exptime=120).download_all()
        for sector in lc_all:
            sectors_all.append(sector.sector)
    elif tic_number in thirty_minute:
        lc_all = eleanor.multi_sectors(tic=tic_number, sectors='all')
        for sector in lc_all:
            sectors_all.append(sector.sector)
    else:
        sectors_all.append('TIC not a JWST standard')
    return sectors_all
 
def lc_1sector(tic_number, sector_number):
    # returns light curve of JWST standard star with TESS data
    # Works for 2min and 30min stars
    # Need two inputs : TIC number and Sector number    
    # The filename variable should be in the local directory where the fits files are stored on your computer.
    
    # example
    #lc = lc_1sector(327587572, 19)
    #lc.scatter()
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    lc = []
    
    if tic_number in two_minute:
        lc.append(lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC',
                                       sector = sector_number, exptime=120).download())
    elif tic_number in thirty_minute:
        filename = ('/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/' + str(tic_number) + 
                    '_sector_' + str(sector_number) + '.fits')
        star=eleanor.Source(tic=tic_number)
        data=eleanor.TargetData(star)
        data.load(directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/', fn=filename)
        lc.append(data.to_lightkurve())
        
    else:
        lc.append('check TIC and sector inputs')
    
    lc2 = lc[0].normalize().remove_nans().remove_outliers()
    
    return lc2

def plot_lightcurve(tic_number, sector_number):
    # Plots normalized lightcurve (in percent); 
    # zero at center with positive and negative values for other data points
    # Works for 2min and 30min stars
    # Need two inputs : TIC number and Sector number    
    # The filename variable should be in the local directory where the fits files are stored on your computer.
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    
    if tic_number in two_minute:
        lc = lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC', 
                                      sector = sector_number, exptime=120).download()
    elif tic_number in thirty_minute:
        filename = ('/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/' + str(tic_number) + 
                    '_sector_' + str(sector_number) + '.fits')
        star=eleanor.Source(tic=tic_number)
        data=eleanor.TargetData(star)
        data.load(directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/', fn=filename)
        lc = data.to_lightkurve()
    else:
        print('TIC not a JWST standard')
        pass
        
    lc2 = lc.remove_nans().remove_outliers().normalize()

    times=lc2.time.value
    fluxes=lc2.flux.value
    
    plt.scatter(times, ((fluxes*100)-100), s=3, c='k', label='TIC {}'.format(tic_number)+
                '\nSector {}'.format(sector_number))
    plt.xlabel('Time [days]', fontsize=15)
    plt.ylabel('Normalized Flux [percent]', fontsize=15)
    plt.legend(loc='lower right', fontsize=11)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.show()
    plt.close()

def upper_limit_percent_variation(tic_number, sector_number):
    
    # Gives the upper limit on percent variation
    # maximum flux - minimum flux
    # This function works for 2min and 30min stars
    # This function should be used only on the variable stars otherwise the percent is meaningless
    # Need two inputs : TIC number and Sector number    
    # The filename variable should be in the local directory where the fits files are stored on your computer.
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    lc = []
    
    if tic_number in two_minute:
        lc.append(lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC',
                                       sector = sector_number, exptime=120).download())
    elif tic_number in thirty_minute:
        filename = ('/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/' + str(tic_number) + 
                    '_sector_' + str(sector_number) + '.fits')
        star=eleanor.Source(tic=tic_number)
        data=eleanor.TargetData(star)
        data.load(directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/', fn=filename)
        lc.append(data.to_lightkurve())
        
    else:
        lc.append('check TIC and sector inputs')
    
    lc2 = lc[0].normalize().remove_nans().remove_outliers()

    times=lc2.time.value
    fluxes=lc2.flux.value
    upperLimit = (np.max(fluxes)-np.min(fluxes))*100
    
    return upperLimit

def flux_variation(tic_number, sector_number):
    # variation of fluxes in light curve in percent units    
    # 99 percentile of fluxes - 1 percentile of fluxes
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    lc = []
    
    if tic_number in two_minute:
        lc.append(lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC',
                                       sector = sector_number, exptime=120).download())
    elif tic_number in thirty_minute:
        filename = ('/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/' + str(tic_number) + 
                    '_sector_' + str(sector_number) + '.fits')
        star=eleanor.Source(tic=tic_number)
        data=eleanor.TargetData(star)
        data.load(directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/', fn=filename)
        lc.append(data.to_lightkurve())
        
    else:
        lc.append('check TIC and sector inputs')
    
    lc2 = lc[0].normalize().remove_nans().remove_outliers()
    fluxes = lc2.flux.value
    variation = (np.percentile(fluxes, 99) - np.percentile(fluxes, 1))*100
    
    return variation

def difference_imaging(tic_number, sector_number, best_period, epoch_time, binsize, tolerance, avg_vmax, diff_vmax):
    # This function works for 2min data
    
    # This function performs a difference image for one sector of one star
    # If the sector does not work try another sector...
    
    # TIC number and sector of the star to perform difference imaging
    # best period.  Should already know the best fit period.  Units = days
    # epoch time shifts the phase folded lightcurve left and right; 
    #            choose a value that plots a max and min like a sine wave so the peak and min can be used 
    # binsize chooses the number of data points to bin; choose a value that makes the phase plot smooth
    # tolerance chooses the amount of cadences (percent of phase curve) 
    #            to average together at the min and max prior to differencing
    # avg_vmax sets the scale for the colorbar of the 'average' plot
    # diff_vmax sets the scale for the colorbar of the 'difference' plot
    # Some good values (but adjust as necessary) that I have used were:
    #             epoch_time=.554, binsize=.008, tolerance=.2, avg_vmax=100000, diff_vmax=1000
    
    tpf = lk.search_targetpixelfile("TIC {}".format(tic_number), sector=sector_number, 
                                    author= 'SPOC', exptime=120).download()
    
    lc = tpf.to_lightcurve(aperture_mask='pipeline')
    lc2 = lc.normalize().remove_nans().remove_outliers()

    best_fit_period = best_period

    folded = lc2.fold(period = best_fit_period, epoch_time=epoch_time)
    
    folded2 = folded.bin(time_bin_size=binsize)
    folded2.plot()

    full_phase_range = folded2.phase[-1].value - folded2.phase[0].value
    tolerance = tolerance * full_phase_range
    min_phase = folded2.time[np.argmin(folded2.flux)].value
    max_phase = folded2.time[np.argmax(folded2.flux)].value

    min_timestamps = folded.time_original[np.where((folded2.time > min_phase - tolerance)
                                                 & (folded2.time < min_phase + tolerance))].value
    max_timestamps = folded.time_original[np.where((folded2.time > max_phase - tolerance)
                                                 & (folded2.time < max_phase + tolerance))].value

    one_quarter_minima = [f for (f, t) in zip(tpf.flux.value, tpf.time.value) if t in min_timestamps]
    one_quarter_maxima = [f for (f, t) in zip(tpf.flux.value, tpf.time.value) if t in max_timestamps]

    avg_image = np.nanmean(tpf.flux.value, axis=0)
    diff_image = np.nanmean(one_quarter_maxima, axis=0) - np.nanmean(one_quarter_minima, axis=0)
    fig, ax = plt.subplots(1,2)
    l = ax[0].imshow(np.flipud(avg_image),cmap = plt.cm.plasma, vmin=0, vmax=avg_vmax) #vmin=0, vmax=100000
    ax[0].set_title('Average Image\nTIC {}'.format(tic_number)+ ' Sector {}'.format(sector_number))
    k= ax[1].imshow(np.flipud(diff_image),cmap = plt.cm.plasma, vmin=0, vmax=diff_vmax) #vmin=0, vmax=5000
    ax[1].set_title('Difference Image\nTIC {}'.format(tic_number)+ ' Sector {}'.format(sector_number))
    fig.set_size_inches((15,6))
    fig.colorbar(l, ax=ax[0])
    fig.colorbar(k, ax=ax[1])
    plt.show()
    plt.close()

def NearbyBrightStars(tic_number):
    # This function works for 2min and 30min stars
    
    # One input: TIC number
    # Lists the nearby bright stars with Tess magnitude < 15 within 200 arcseconds
    # First object in list is the input star...
    # Object info listed: TIC, Tmag, Jmag, Teff, logg, object type, distance [arcseconds] from input TIC star
    
    radSearch = .056 #  <-- radius in degrees = 200 arcseconds
    
    catalogData = Catalogs.query_object('TIC {}'.format(tic_number), radius = radSearch, catalog = "TIC")

    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']

    # Create a list of nearby bright stars (tess magnitude less than 14) from the rest of the data for later.
    bright = catalogData['Tmag'] < 15

    # Make it a list of Ra, Dec pairs of the bright ones. This is now a list of nearby bright stars.
    #nearbyStars_coordinates = list( map( lambda x,y:[x,y], catalogData[bright]['ra'], catalogData[bright]['dec'] ) )
    nearbyStars = catalogData[bright]['ID', 'Tmag', 'Jmag', 'Teff', 'logg', 'objType', 'dstArcSec']

    return nearbyStars
 
def saving_eleanor_data_to_fits_file_locally(tic_number):
    
    # Two minute data is not avalailable for all of the stars.
    # For those without 2minute data we use full frame images which have thirty minute data.
    # The thirty minute data has many systematic errors so we use eleanor to obtain the fits files of corrected data.

    # For stars with 30 minute data:
    # The length of the thirty_minute array is 14 (13 is last index), ie 14 stars to get FFIs for.
    # This loop saves each fits file on disk for each star and each available sector.
    # There are 112 fits files total.

    # This loop saves all sectors of data for one star.
    # The filename changes automatically.  Saved as TICnumber_sector_Sectornumber.fits, 
    # where TIC number and Sector number would be changing.

    # Change the directory to where you want the fits files saved locally.

    star=tic_number   
    star1 = eleanor.multi_sectors(tic=star, sectors='all')

    i = 0

    while i < len(star1):

        starsector = star1[i]
        data = eleanor.TargetData(starsector, height=15, width=15, bkg_size=31, do_psf=True, do_pca=True, 
                                  regressors='corner')
        data.save(output_fn='{}'.format(star) + '_sector_{}.fits'.format(starsector.sector), 
                  directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz')

        i+=1 
def see_background_stars_in_tpf(tic_number, sector_number):
    # Plots target pixel file with all of the background stars and their stats as well,
    # so that we can see their distance/size/magnitude in comparison to our target object
    # This function works for 2min data
    lc = lk.search_targetpixelfile('TIC {}'.format(tic_number), 
                                   sector=sector_number, author= 'SPOC', exptime=120).download()
    lc.interact_sky() 

####################################################
#import k2flix
#Converts a Target Pixel File (TPF) from NASA's Kepler/K2/TESS spacecraft into
#an animated gif or MPEG-4 movie for human inspection.
#from console:
#   to make a gif:
#        k2flix tpf-file.fits.gz       
#   to make a mpeg4 movie:
#        k2flix --o movie.mp4 tpf-file.fits.gz
####################################################

def two_or_thirty_or_none(tic_number):
    
    # tells you whether the star has 2min, 30min, or no data
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    none = [352817378,135656809,397558558,420814525,54837036,315627636,80313923,75586606]
    statement = []
    if tic_number in two_minute:
        statement.append('two minute')
    elif tic_number in thirty_minute:
        statement.append('thirty minute')
    elif tic_number in none:
        statement.append('Neither 2 nor 30 min data')
    else:
        statement.append('TIC input was not a JWST standard')
    return statement[0]

def info_on_single_star(tic_number):
    # This function works for the 2min and 30min stars
    # It gives back information about the star with input TIC
    # One input for function: TIC number
    # Object info listed: TIC, ra, dec, Tmag, Vmag, Kmag, Teff, logg
    
    radSearch = .056 #  <-- radius in degrees = 200 arcseconds
    
    catalogData = Catalogs.query_object('TIC {}'.format(tic_number), radius = radSearch, catalog = "TIC")

    info = catalogData[0]['ID', 'ra', 'dec', 'Tmag', 'Vmag', 'Kmag', 'Teff', 'logg']

    return info

def info_on_all_stars(name_of_csv):
    # This function will create a csv in excel format with info about ALL of the JWST standards
    # Input is the name of the csv file; variable must be a string with .csv at the end
        # example: 'statistics.csv'
    # Info returned: TIC, ra, dec, Tess magnitude, V mag, K mag, effective temperature, and log G
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    none = [352817378,135656809,397558558,420814525,54837036,315627636,80313923,75586606]
    variable_list = [383553764,166698220,41232189,198456033,441120034,32869782,219820925]
    full_list_of_JWST_standards = two_minute+thirty_minute+none
    
    with open(name_of_csv, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', dialect='excel')
        spamwriter.writerow(['ID', 'ra', 'dec', 'Tmag', 'Vmag', 'Kmag', 'Teff', 'logg']) 

        radSearch = .056 #  <-- radius in degrees = 200 arcseconds

        for star in full_list_of_JWST_standards:
            radSearch = .056 #  <-- radius in degrees = 200 arcseconds

            catalogData = Catalogs.query_object('TIC {}'.format(star), radius = radSearch, catalog = "TIC")

            info = catalogData[0]['ID', 'ra', 'dec', 'Tmag', 'Vmag', 'Kmag', 'Teff', 'logg']


            a = info['ID']
            b = info['ra']
            c = info['dec']
            d = info['Tmag']
            e = info['Vmag']
            f = info['Kmag']
            g = info['Teff']
            h = info['logg']

            spamwriter.writerow([a] + [b] + [c] + [d] + [e] + [f] + [g]+[h])
        
def plot_periodogram(tic_number, sector_number, max_freq):
    
    # Need 3 inputs: 
    #       -TIC ID 
    #       -Sector Number
    #       -Maximum frequency to plot up to on the x-axis [in uHz]
    
    # Plots the lombscargle periodogram 
    # Y axis = amplitude [percent ie parts per hundred]
    # bottom X axis = frequency [microhertz ie. uHz]
    # upper X axis = period [hours]
    
    # Function works for 2min and 30min stars.
    # For 30 min data, the filename variable should be in the local directory 
    #     where the fits files are stored on your computer.
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    
    if tic_number in two_minute:
        lc = lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC', 
                                      sector = sector_number, exptime=120).download()
    elif tic_number in thirty_minute:
        filename = ('/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/' + str(tic_number) + 
                    '_sector_' + str(sector_number) + '.fits')
        star=eleanor.Source(tic=tic_number)
        data=eleanor.TargetData(star)
        data.load(directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/', fn=filename)
        lc = data.to_lightkurve()
    else:
        print('TIC not a JWST standard')
        pass
        
    lc2 = lc.remove_nans().remove_outliers().normalize()
    
    pgram=lc2.to_periodogram(method='lombscargle', normalization='amplitude')
    
    N=len(pgram.frequency)
    
    sigma_rms = np.nanstd(lc2.flux)
    # mean noise level in amplitude spectrum
    sigma_amp = np.sqrt(np.pi/N)*sigma_rms
    
    freq=pgram.frequency
    amp=pgram.power
    
    microhertz=11.574074074*freq
    hertz = freq*0.000011574074074
    per_hour = freq*0.041666666667
    
    fig, ax1 = plt.subplots()

    ax1.plot(microhertz, amp*100, color='k', label='TIC {}'.format(tic_number)+
             '\nSector {}'.format(sector_number))
    plt.xlim(0, max_freq)
    plt.xticks(fontsize=15)
    plt.yticks([0,.1,.2,.3,.4], fontsize=15)
    plt.xlabel('Frequency [uHz]', fontsize=15)
    plt.ylabel('Amplitude [percent]', fontsize=15)
    plt.hlines(5*sigma_amp*100, 0, max_freq, colors='r', linestyles='dashed', label='5 x Mean\nNoise Level')
    plt.legend(loc='upper right', fontsize=10)

    ax2 = ax1.twiny()
    plt.xlabel('Period [hours]', fontsize=15)
    ax2.set_xticklabels([round(24*11.574074074/x, 1) for x in ax1.get_xticks()])
    plt.xticks(fontsize=15)
    plt.show()
    plt.close()

def fold(tic_number, sector_number):
    # Plots folded, normalized lightcurve (in percent); 
    # zero at center with positive and negative values for other data points
    # Works for 2min and 30min stars
    # Need two inputs : TIC number and Sector number    
    # The filename variable should be in the local directory where the fits files are stored on your computer.
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    
    if tic_number in two_minute:
        lc = lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC', 
                                      sector = sector_number, exptime=120).download()
    elif tic_number in thirty_minute:
        filename = ('/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/' + str(tic_number) + 
                    '_sector_' + str(sector_number) + '.fits')
        star=eleanor.Source(tic=tic_number)
        data=eleanor.TargetData(star)
        data.load(directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/', fn=filename)
        lc = data.to_lightkurve()
    else:
        print('TIC not a JWST standard')
        pass
        
    lc2 = lc.remove_nans().remove_outliers().normalize()
    
    pgram=lc2.to_periodogram(method='lombscargle', normalization='amplitude')
    best_fit_period = pgram.period_at_max_power
    
    lc3 = lc2.fold(period=best_fit_period, normalize_phase=True)
    
    times=lc3.time.value
    fluxes=lc3.flux.value
    
    plt.scatter(times, (fluxes*100-100), s=3, c='k', label='TIC {}'.format(tic_number)+
                '\nSector {}'.format(sector_number))
    plt.xlabel('Phase', fontsize=15)
    plt.ylabel('Normalized Flux [percent]', fontsize=15)
    plt.legend(loc='lower right', fontsize=11)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.show()
    plt.close()

def crowd_sap_1star_1sector(tic_number, sector_number):
    # returns crowd sap value of one sector of one star
    # Function works for 2min data
    # example:
    lc = lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC', 
                              sector = sector_number, exptime=120).download()
    return lc.crowdsap

def crowd_sap_all_2min(array):
    
    # Function works for 2min data
    # saves all crowd saps for all sectors of all 2min stars to csv in excel format
    # example:
    #     crowd_sap(two_minute)
    
    with open('crwdsaps.csv', 'w', newline='') as csvfile:

        spamwriter = csv.writer(csvfile, delimiter=',', dialect='excel')
        spamwriter.writerow(['TIC', 'sector', 'crowd sap']) 

        for tic_number in array:
            lc = lk.search_lightcurve('TIC {}'.format(tic_number), 
                                      author='SPOC', exptime=120).download_all()
            z=0
            while z < len(lc):
                a = tic_number
                b = lc[z].sector
                c = lc[z].crowdsap
                
                spamwriter.writerow([a] + [b] + [c]) 
                
                z+=1                

def bootstrap(tic_number, sector_number):
    
    # Getting the noise level for 1 sector of 1 star; function returns 5 x avg noise level in percent!
        
    # one in a thousand false alarm probability
    # output is 5 times the noise level reported from previous equation 
        # previous equation -> sigma_amp = np.sqrt(np.pi/N)*sigma_rms, 
        # where sigma_amp = 1 * noise level, sigma_rms = stdev of fluxes, and N = number of data points

    # the fct does the following:
        # takes all of the times and fluxes, 
        # shuffle the fluxes,
        # takes a periodogram, 
        # measures the highest peak, 
        # keeps track of it, 
        # performs this 10K times
        
    # ie does 10,000 shuffles
    # report "99 % of the time never see a peak higher than _____ due to noise alone"
    
    # since the fluxes are shuffled the time series has no periodicity
    # therefore the max amplitude of the periodogram is due to noise alone
    
    # period range for the periodogram is 2 hours (.0833 days) to 10 days
    # max_amp_list = list of the highest amps for each of the 10K shuffles
    # np.percentile(max_amp_list, 99.9)
    
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    
    if tic_number in two_minute:
        lc = lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC', 
                                      sector = sector_number, exptime=120).download()
    elif tic_number in thirty_minute:
        filename = ('/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/' + str(tic_number) + 
                    '_sector_' + str(sector_number) + '.fits')
        star=eleanor.Source(tic=tic_number)
        data=eleanor.TargetData(star)
        data.load(directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/', fn=filename)
        lc = data.to_lightkurve()
    else:
        print('TIC not a JWST standard')
        pass
    
    #########################################    
    def shuffle(aList):
        return random.shuffle(aList)
    ###########################################
    
    lc2 = lc.remove_nans().remove_outliers().normalize()
    times=lc2.time.value
    fluxes=lc2.flux.value
    
    max_amp_list = []
    N = 10000 # number of shuffles
    n = 0
    while n < N:

        shuffle(fluxes)

        lc3 = lk.LightCurve(time = times, flux = fluxes)

        pgram=lc3.to_periodogram(method='lombscargle', normalization='amplitude', 
                                 minimum_period = .0833, maximum_period = 10)
        max_amp = pgram.max_power*100 # times 100 to turn into percent
        max_amp_list.append(max_amp)

        n+=1
        
    # use 99.9 percentile because 10K shuffles and 1/1000 probability of amp = 5*(avg noise level)
    five_times_noise_level = np.percentile(max_amp_list, 99.9)
    return five_times_noise_level

def stitch(tic_number):
    
    # makes a light curve of ALL available sectors for one star
    # works for 2min and 30 min data
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    
    times = []
    fluxes = []
    
    if tic_number in two_minute:
        lc_all = lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC', exptime=120).download_all()
        for lc in lc_all:
            lc2 = lc.normalize().remove_nans().remove_outliers()
            times.extend(lc2.time.value)
            fluxes.extend(lc2.flux.value)
    elif tic_number in thirty_minute:
        lc_all = eleanor.multi_sectors(tic=tic_number, sectors='all')
        for sector in lc_all:
            sector_number = sector.sector
            filename = ('/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/' + str(tic_number) + 
                        '_sector_' + str(sector_number) + '.fits')
            star=eleanor.Source(tic=tic_number)
            data=eleanor.TargetData(star)
            data.load(directory='/Users/mkunz/Tess_JWST_calibration_stars/thirty_minz/', fn=filename)
            lc = data.to_lightkurve()
            lc2 = lc.normalize().remove_nans().remove_outliers()
            times.extend(lc2.time.value)
            fluxes.extend(lc2.flux.value)
    else:
        pass
    lc3 = lk.LightCurve(time=times, flux=fluxes)
    lc3.scatter(label='TIC {}'.format(tic_number))

def one_noise_level(N, fluxes):
    
    # Funtion returns 1*average noise level (in percent!) in time series data
    # This is an alternative to using the bootstrap method
    # N = len(lc.time) ; ie number of data points
    # sigma_amp = np.sqrt(np.pi/N)*sigma_rms, 
    # where sigma_amp = 1 * noise level, sigma_rms = stdev of fluxes, and N = number of data points
    
    # Example of using this function:
    # lc = lc_1sector(41232189, 3)
    # N = len(lc.time)
    # fluxes = lc.flux.value
    # noise = one_noise_level(N, fluxes)
    # print(noise)
    
    sigma_rms = np.std(fluxes) # root mean square is standard deviation of flux values
    sigma_amp = np.sqrt(np.pi/N)*sigma_rms
    return sigma_amp*100 # multiply by 100 to be in percent units

def bootstrap_all_standards(name_of_csv):
    # bootstraps all JWST standards and saves the data to csv in excel format
    # info reported is TIC, sector, and 5*NoiseLevel
    # Input is name of csv file
    # Example:
    # name_of_csv = 'noise_levels.csv'
    # bootstrap_all_standards(name_of_csv)
    
    practice_list = [147921014]
    two_minute = [327587572,247923021,149505899,147921014,471015233,383553764,166698220,
                  41232189,298165335,198456033,181240911,165370459,219752116,27533327,
                  39464221,441120034,140282069,32869782,365653206,229980646]
    thirty_minute = [219820925,229945862,440765193,8591766,417544924,144599609,207440438,
                     219094190,233095291,219114641,233067231,233075513,219897252,233205654]
    none = [352817378,135656809,397558558,420814525,54837036,315627636,80313923,75586606]
    variable_list = [383553764,166698220,41232189,198456033,441120034,32869782,219820925]
    full_list_of_JWST_standards = two_minute+thirty_minute+none
    
    with open(name_of_csv, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', dialect='excel')
        spamwriter.writerow(['TIC', 'Sector', 'Bootstrap [5*NoiseLevel]']) 
        
        # practice_list is the array of tic numbers 
        #      that the function will loop through to report data on
        # can change practice_list to be two_minute, thirty_minute, or full_list_of_JWST_standards
        
        for tic_number in practice_list:
            if tic_number in two_minute:
                lc_all = lk.search_lightcurve('TIC {}'.format(tic_number), author='SPOC', 
                                              exptime=120).download_all()
                for lc in lc_all:
                    sector_number = lc.sector
                    a = tic_number
                    b = sector_number
                    c = bootstrap(tic_number, sector_number)
                    spamwriter.writerow([a] + [b] + [c])
            elif tic_number in thirty_minute:
                lc_all = eleanor.multi_sectors(tic=tic_number, sectors='all')
                for lc in lc_all:
                    sector_number = lc.sector
                    a = tic_number
                    b = sector_number
                    c = bootstrap(tic_number, sector_number)
                    spamwriter.writerow([a] + [b] + [c])
            else:
                pass

