# -*- coding: utf-8 -*-
"""
stat_plot_csv.py

This file produces either a .csv file of statistics for the output files of a particular source
or a series of plots (saved as PDFs) for each visit for the source. To select which to perform, 
the user executes the code through the command line using 'python stat_plot.py x', where 'x'
is either 'stats' for the statistics file, 'plots' for the PDF plots, 'movies' for .fits file 
data cubes, or 'fchart' for finding charts, or a .fits image of the source. 

This file utilizes existing files for its output, namely the time, find and CSV files for 
each source and its respective visits. The time files are exact output from using the gFind
module to search for sources; they reveal the number of exposures in a certain band as well as
the time ranges in which these exposures occur. The find files are files that are created 
before launching gFind and contain information relevant to the source, such as RA, Dec, 
and FUV & NUV mags. The CSV files comprise photometry for a source in a certain band during 
a specific visit, along with other parameters. This file incorporates all of these types of
files to concisely spit out the type of file the user desires. 

Alexander de la Vega, Aggie Quesada (JHU)
Latest version: 2 May 2017
"""

from __future__ import division # python 2.x usually has issues dividing
import numpy as np # for computation / array formatting
import os # move through directories
import sys # read in command line arguments
from astropy.io import ascii


try:
    os.environ['GPHOTON_HOME']
except KeyError:
    print "gPhoton environment variables not found. You need to either:\n 1. Source ~/path/to/gphoton/home/gphoton_var.sh in this session and try again, or\n 2. Add 'Source ~/path/to/gphoton/home/gphoton_var.sh' to your ~/.bash_profile"
    sys.exit(0)

def is_float(x): # check whether values from CSVs are floats or not
    try:
        float(x)
        return True
    except ValueError:
        return False

time_path = os.environ['GPHOTON_TIMES'] # path to store time ranges for each target
if not os.path.exists(time_path): # check whether times directory exists
    os.makedirs(time_path) # if not, create directory

files = os.listdir(time_path) # list of all files in times directory (i.e. each object of interest)
num_file = len(files) # number of files in times directory

csv_path = os.environ['GPHOTON_OUTPUT']
if not os.path.exists(csv_path): # check whether output directory exists
    os.makedirs(csv_path) # if not, create directory

files_csv = os.listdir(csv_path) # list of CSVs in CSV directory
num_csv = len(files_csv) # number of CSVs

out_path = os.environ['GPHOTON_PLOTS'] # output path
if not os.path.exists(out_path): # check whether output directory exists
    os.makedirs(out_path) # if not, create directory

find_path = os.environ['GPHOTON_FIND'] # path to store gFind output for each target
if not os.path.exists(find_path): # check whether gFind output directory exists
    os.makedirs(find_path) # if not, create directory

if len(sys.argv) > 2: # if input is given as .txt file
    if sys.argv[2].endswith('.txt'):
        inputf = open(sys.argv[2], 'r')
        input_cont = inputf.readlines()
        num_file = len(input_cont)
        file_list = [n.replace('\n', '') for n in input_cont]
        inputf.close()

aperture = str(os.environ['GPHOTON_APERTURE']) # aperture radius in arcsec
time_step = str(os.environ['GPHOTON_TIME_BIN']) # time step in sec

aper = str(os.environ['GPHOTON_APERTURE']) # vector of aperture radii strings in arcsec; used to construct file names for output files
ts = str(os.environ['GPHOTON_TIME_BIN']) # time step in sec - vector of strings for file name construction mentioned above
band = ['FUV', 'NUV'] # vector of photometric bands - for aforementioned file name construction

# if user passes 'stats' to the command line - produce stats output .txt file  
if sys.argv[1] == 'stats':
    from scipy.signal import argrelextrema
    from astropy.convolution import convolve_fft, Box1DKernel

    kernel = Box1DKernel(3)
    char_freq = 0.00882353

    good_flag = '0.0'

    stats_file = open(out_path + 'fft_stats.csv', 'w')
    stats_file.write('objid,band,visit,duration,min_mag,max_mag,mean_mag_bgsub,max_mag_var,max_mag_var_sigma,flag0_max_mag_var,flag0_max_mag_var_sigma,max_counts_var_ratio,max_counts_var_err_ratio,' +                 'mean_det_rad,err_det_rad,num_harmonics,num_other_peaks,flags,frac_bins_flags,tot_frac_bins_flags\n')

    small_exptime_file = open(out_path + 'low_exposure_time_files.txt', 'w')
    miss_val_file = open(out_path + 'missing_value_files.txt', 'w')

    file_list = os.listdir(csv_path)

    for f in file_list:
        if f.endswith('.csv'):
            objid = f.split('_')[0]
            band = f.split('_')[7]
            visit = f.replace('.csv', '').split('_')[8][4:]
            csv_f = ascii.read(csv_path + f)
            
            duration = csv_f['t0'][len(csv_f) - 1] - csv_f['t0'][0]
            
            csv_f['t_mean'].fill_value = -99
            csv_f['detrad'].fill_value = -99
            
            try:
                if -99 in csv_f['t_mean'].filled() or -99 in csv_f['detrad'].filled():
                    print f
                    miss_val_file.write(f + '\n')
                    continue
            except AttributeError:
                pass
            
            good_exp = np.where(csv_f['exptime'] >= 0.75 * time_step)[0]
            if len(good_exp) == 0:
                #print f
                small_exptime_file.write(f + '\n')
                continue
            
            min_mag = np.min(csv_f['mag_bgsub'][good_exp])
            max_mag = np.max(csv_f['mag_bgsub'][good_exp])
            
            mean_mag_bgsub = np.mean(csv_f['mag_bgsub'][good_exp])
            mean_mag_bgsub_err = np.mean(csv_f['mag_bgsub_err_2'][good_exp])
            
            mean_det_rad = np.mean(csv_f['detrad'][good_exp])
            err_det_rad = np.std(csv_f['detrad'][good_exp])
            
            max_mag_var = np.max(abs(csv_f['mag_bgsub'][good_exp] - mean_mag_bgsub))
            max_mag_index = np.where(abs(csv_f['mag_bgsub'][good_exp] - mean_mag_bgsub) == max_mag_var)
            max_mag_var_sigma = max_mag_var / csv_f['mag_bgsub_err_2'][max_mag_index]
            
            mean_counts = np.mean(csv_f['flat_counts'][good_exp])
            max_counts_var = np.max(abs(csv_f['flat_counts'][good_exp] - mean_counts))
            max_counts_index = np.where(abs(csv_f['flat_counts'][good_exp] - mean_counts) == max_counts_var)[0]
            max_counts_var_ratio = max_counts_var / mean_counts
            max_counts_var_err_ratio = csv_f['cps_err'][good_exp][max_counts_index] *            csv_f['exptime'][good_exp][max_counts_index] / mean_counts
            
            data_fft = scpfft.fft(csv_f['mag_bgsub'])
            n = len(csv_f['mag_bgsub'])
            freqs = np.linspace(1, .5 * n - 1, .5 * n - 1) * .2 / n
            fft_half = abs(data_fft[1:int(.5 * n)])
            
            fft_conv = convolve_fft(fft_half, kernel)
            maxima = argrelextrema(fft_conv, np.greater)[0]
            
            freq_ratios = np.array([freqs[i] for i in maxima]) / (char_freq)
            
            max_num_harmonics = int(((.5 * n - 1) * .2 / n) / char_freq)
            num_harmonics = funcs.find_num_harmonic(freq_ratios)
            num_pot_var = len(maxima) - num_harmonics
            
            if num_harmonics > max_num_harmonics:
                num_harmonics = max_num_harmonics
                
            flags = str(funcs.find_unique_elements(csv_f['flags'][good_exp])).replace(',', '')
            
            if check_flag(flags, good_flag):
                good_flag0_criteria = (csv_f['flags'] == float(good_flag)) #& (csv_f['exptime'] >= 0.75 * time_step) 
                
                flag0_mean_mag_bgsub = np.mean(csv_f['mag_bgsub'][good_flag0_criteria])
                flag0_mean_mag_bgsub_err = np.mean(csv_f['mag_bgsub_err_2'][good_flag0_criteria])
                
                flag0_max_mag_var = np.max(abs(csv_f['mag_bgsub'][good_flag0_criteria] - flag0_mean_mag_bgsub))
                flag0_max_mag_var_sigma = flag0_max_mag_var / flag0_mean_mag_bgsub_err
            else:
                flag0_max_mag_var = -99.0
                flag0_max_mag_var_sigma = -99.0
                
            tot_frac_bins_flags = '%.2f ' % (len(csv_f['flags'][good_exp][csv_f['flags'][good_exp] > 0]) / len(csv_f[good_exp]))
            frac_bins_flags = [len(csv_f[csv_f['flags'] == i]) for i in funcs.find_unique_elements(csv_f['flags'][good_exp])]
            frac_bins_flags = np.asarray(frac_bins_flags) / len(csv_f[good_exp])
            frac_bins_flags = ['%.2f' % i for i in frac_bins_flags]
            frac_bins_flags = str(frac_bins_flags).replace(',', '').replace("'", "")
            
            stats_file.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19}\n'.format(objid,
                band,visit,str(duration),str(min_mag),str(max_mag),str(mean_mag_bgsub),str(max_mag_var),str(max_mag_var_sigma),\
                str(flag0_max_mag_var),str(flag0_max_mag_var_sigma),str(max_counts_var_ratio),\
                str(max_counts_var_err_ratio[0]),str(mean_det_rad),str(err_det_rad),str(num_harmonics),\
                str(num_pot_var),flags,frac_bins_flags,tot_frac_bins_flags))

    stats_file.close()
    small_exptime_file.close()
    miss_val_file.close()


# if the user passes 'plots' to the command line
elif sys.argv[1] == 'plots':
    import matplotlib.pylab as plt # plotting
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter # for tick marks
    from matplotlib.backends.backend_pdf import PdfPages
    import datetime
    from pylab import rcParams
    import glob

    py_path = '/Users/Aggie/Documents/Research/CSVpy/CSVpy/' #'./CSVpy/' #path storing the csv.py files
    light_path = '/Users/Aggie/Documents/Research/test/' #output path desired for the plots
    csv_path = '/Users/Aggie/Documents/Research/Data/out/'  #path storing the csvs 
    
    if not os.path.exists(py_path):
        print 'Please check input for CSVpy path'
        sys.exit(0)

    file_list_csvpy = os.listdir(py_path)
    
    def plotfile(csvpy_f, aper): #pick aper 15 or 25
    
        # From files need magnitudes, errors, exposure time. Open the files and store everything   
        csvpy_path = py_path + csvpy_f
        csvpy_file = open(csvpy_path, 'r')
        content = csvpy_file.readlines()# Reads all the lines of the file and returns them into a list
        csvpy_file.close()

        #Find where the line telling the exposure time is
        for i in range(len(content)):
            if len(content[i]) > 1: # If not an empty line
                if content[i].split()[0] == '#FUV': # This only shows up once in the file
                    f_line = i 
                if content[i].split()[0] == '#NUV':
                    n_line = i
                    break  # Cut the loop. Fuv guaranteed to be before nuv (otherwise it would break)

        # Get the exposure times (number of visits)
        f_exp, n_exp = int(content[f_line].split()[3]), int(content[n_line].split()[3]) 


        #Store the begin and end times of the exposures into arrays
        begin_times_f, end_times_f = np.empty(f_exp), np.empty(f_exp) #np.empty creates an empty array of that size
        begin_times_n, end_times_n = np.empty(n_exp), np.empty(n_exp)

        for i in range(f_line + 1, f_line + f_exp + 1): #+1 b/c range() is noninclusive
            #begin_times is a array that we want to map from the lines in the file. do [i-(f_line +1))] in order to map correctly since 
            #f_line isn't necessarily 0 to start with and the indicies in the array actually do start with 0 
            begin_times_f[i - f_line - 1] = float(content[i].replace('[','').replace(']','').replace(',','').split()[1]) 
            end_times_f[i - f_line - 1] = float(content[i].replace('[','').replace(']','').split()[2])

        for j in range(n_line + 1, n_line + n_exp + 1):
            begin_times_n[j - n_line - 1] = float(content[j].replace('[','').replace(']','').replace(',','').split()[1])
            end_times_n[j - n_line - 1] = float(content[j].replace('[','').replace(']','').split()[2])

        #Find where the visits of NUV and FUV match up
        time_matches_f = np.asarray([i for i in range(1, f_exp + 1)]) # i for i in range creates a list. np.asarray turns a list into an array

        if n_exp >= f_exp:
            time_matches_n = np.zeros(n_exp)
            count = 0
            for k in range(n_exp):
                if begin_times_n[k] in begin_times_f:
                    count += 1
                    time_matches_n[k] = count #marks as the kth element as the countith match
                    

        #If something in time_matches_n = 0, no NUV during that visit
        #Time_matches_f doesn't get used elsewhere. Might have just been for diagnostics

        #Find min and max mag limits to determine the axes in the light plots below 


        #Now let's plot each visit
        #the below is for cosmetics
        params = {'text.usetex': False, 'mathtext.fontset': 'stixsans'} #turn off latex labels, change math font 
        plt.rcParams.update(params)  #rc file does permanent configuration change on plt

        #In order to determine y-axes in the light curves: Find min and max magnitude limits
        file_list = glob.glob(csv_path + '*.csv')
        min_mag_limit = 20
        max_mag_limit = 5
        minlims = []
        maxlims = []
        for f in file_list:
            filename = os.path.split(f)[1]
            if filename.startswith(csvpy_f.split('_')[0]):
                csv_f = ascii.read(f)
                f_min_mag = np.min(csv_f['mag_bgsub'] - csv_f['mag_bgsub_err_1'])
                f_max_mag = np.max(csv_f['mag_bgsub'] + csv_f['mag_bgsub_err_2'])
                if f_min_mag < min_mag_limit:
                    minlims.append(f_min_mag)
                if f_max_mag > max_mag_limit:
                    maxlims.append(f_max_mag)
        min_mag_limit = np.min(minlims)
        max_mag_limit = np.max(maxlims)
        
        #so that the y axis isn't too small:
        if max_mag_limit - min_mag_limit <= 1:
            max_mag_limit += 0.5
            min_mag_limit -= 0.5

        # 9.0 comes from available real estate on page.
        # 1.0 adjustment factor comes from experimentation.
        numPlotsPerPage = int(9.0 / (max_mag_limit - min_mag_limit) * 1.0)
        if numPlotsPerPage == 0: numPlotsPerPage = 1 # If there's a crazy-big range, int() can yield 0. Fix that here!


        #Because some data only has NUV observations, we need to check if there are visits that match up. 
        #If only NUV, plot only NUV. If both FUV and NUV, plot both
        if n_exp >= f_exp:

            #Create subplots so I can save all the visits for each source under one image. 
            with PdfPages(light_path + '%s_%d.pdf' % (csvpy_f.split("_")[0], aper)) as pdf:
                fig, subplots = plt.subplots(numPlotsPerPage, figsize=(10,12))
                if numPlotsPerPage == 1: subplots = [subplots] # Before this, subplots is not an array. Fix to avoid indexing errors down the line
                n = 0 # Number of plots on page right now

                for i in range(n_exp):

                    if time_matches_n[i] == 0:

                        #Read the appropriate NUV file
                        fname = csvpy_f.split('_')[0] + '_aper_' + str(aper) + '_as_timestep_5_sec_NUV_time' + str(i + 1) + '.csv'
                        csv_f = ascii.read(csv_path + fname)


                        #Find location on the detector
                        xloc = np.mean(csv_f['detxs']) * 0.1 #units were in 0.1 arcmin
                        yloc = np.mean(csv_f['detys']) * 0.1 

                        #Find detector radius
                        rad = np.mean(csv_f['detrad']) * 0.1 

                        #Find timestamp
                        value = datetime.datetime.fromtimestamp(csv_f["t0"][0] + 315964800) #offset between galex time step and unix
                        timestamp = value.strftime('%Y-%m-%d %H:%M:%S') 


                        #X-Axes: Find the correct time to put on the x-axis
                        time_range = csv_f['t0'][len(csv_f) - 1] - csv_f['t0'][0]
                        time_buffer = 0.05 * time_range


                        #Distinguish between good and bad criteria
                        good_criteria = (csv_f['exptime'] >= .75 * 5)
                        bad_criteria = (csv_f['exptime'] < .75 * 5)


                        #Set minor tick marks based on length of visit and set plot size
                        if time_range <= 200: # if length <= 200 sec
                            minorLocator = MultipleLocator(5) # tick marks every 5 sec
                        else:
                            minorLocator = MultipleLocator(20) # tick marks every 20 sec

                        subplots[n].xaxis.set_minor_locator(minorLocator) # plot with minor tick marks  
                        subplots[n].tick_params(axis='both', which='major', labelsize=8)

                        #Plot only NUV plot
                        #2 y-errors bc one for above and one for below
                        subplots[n].errorbar(csv_f['t0'][good_criteria] - csv_f['t0'][0] + time_buffer, csv_f['mag_bgsub'][good_criteria],                                        yerr=[csv_f['mag_bgsub_err_2'][good_criteria], csv_f['mag_bgsub_err_1'][good_criteria]],                                         linestyle='None',marker='o', markersize=3, color='r')

                        #Bad criteria with markerfacecolor='None' makes the dot have no filling
                        subplots[n].errorbar(csv_f['t0'][bad_criteria] - csv_f['t0'][0] + time_buffer, csv_f['mag_bgsub'][bad_criteria],                                        yerr=[csv_f['mag_bgsub_err_2'][bad_criteria], csv_f['mag_bgsub_err_1'][bad_criteria]],                                        linestyle='None', marker='o', markersize=3, markerfacecolor='None', color='r')


                        #Axes
                        subplots[n].set_ylim(max_mag_limit + 0.05, min_mag_limit - 0.05)
                        subplots[n].set_xlim(0, time_buffer * 22)
                        subplots[n].set_ylabel('AB mag') #AB magnitude = astronomical magnitude system--since magnitudes are relative this determines 0 as something set (eg: Vega star brightness)
                        subplots[n].set_xlabel('Time - ' + str('%.f' % (csv_f['t0'][0] - time_buffer)) + ' (s)')

                        #Title
                        param = fname.replace('.csv', '').replace('time','').split('_')
                        subplots[n].set_title(param[0] + ' ' + param[7] + ' Visit ' + param[8] + ' - aper ' + param[2] + ' as, timestep ' + param[5] + (' s\nDetloc x/y %.2f/%.2f arcmin, Detrad %.2f arcmin' % (xloc, yloc, rad)) + ', Time ' + timestamp)

                    else:


                        #Read the appropriate NUV and FUV files
                        fname_n = csvpy_f.split('_')[0] + '_aper_' + str(aper) + '_as_timestep_5_sec_NUV_time' + str(i + 1) + '.csv'
                        fname_f = csvpy_f.split('_')[0] + '_aper_' + str(aper) + '_as_timestep_5_sec_FUV_time' + str(int(time_matches_n[i])) + '.csv'
                        csv_f_n = ascii.read(csv_path + fname_n)
                        csv_f_f = ascii.read(csv_path + fname_f)


                        #Find location on the detector
                        xloc = np.mean(csv_f['detxs']) * 0.1 #units were in 0.1 arcmin
                        yloc = np.mean(csv_f['detys']) * 0.1 

                        #Find detector radius
                        rad = np.mean(csv_f['detrad']) * 0.1 

                        #Find timestamp
                        value = datetime.datetime.fromtimestamp(csv_f["t0"][0] + 315964800) #offset between galex time step and unix
                        timestamp = value.strftime('%Y-%m-%d %H:%M:%S') 


                        #X-Axes: Find the correct time to put on the x-axis
                        if len(csv_f_f) >= len(csv_f_n):
                            time_range = csv_f_f['t0'][len(csv_f_f) - 1] - csv_f_f['t0'][0]
                        else:
                            time_range = csv_f_n['t0'][len(csv_f_n) - 1] - csv_f_n['t0'][0]

                        time_buffer = 0.05 * time_range


                        #Distinguish between good and bad criteria
                        good_criteria_f = (csv_f_f['exptime'] >= .75 * 5)
                        bad_criteria_f = (csv_f_f['exptime'] < .75 * 5)

                        good_criteria_n = (csv_f_n['exptime'] >= .75 * 5)
                        bad_criteria_n = (csv_f_n['exptime'] < .75 * 5)


                        #Set minor tick marks based on length of visit and set plot size
                        if time_range <= 200: 
                            minorLocator = MultipleLocator(5) 
                        else:
                            minorLocator = MultipleLocator(20) 


                        subplots[n].xaxis.set_minor_locator(minorLocator) # plot with minor tick marks 
                        subplots[n].tick_params(axis='both', which='major', labelsize=8)


                        #fuv plot
                        subplots[n].errorbar(csv_f_f['t0'][good_criteria_f] - csv_f_f['t0'][0] + time_buffer, csv_f_f['mag_bgsub'][good_criteria_f],                                        yerr=[csv_f_f['mag_bgsub_err_2'][good_criteria_f], csv_f_f['mag_bgsub_err_1'][good_criteria_f]],                                         linestyle='None', marker='o', markersize=3)

                        subplots[n].errorbar(csv_f_f['t0'][bad_criteria_f] - csv_f_f['t0'][0] + time_buffer, csv_f_f['mag_bgsub'][bad_criteria_f],                                        yerr=[csv_f_f['mag_bgsub_err_2'][bad_criteria_f], csv_f_f['mag_bgsub_err_1'][bad_criteria_f]],                                         linestyle='None', marker='o', markersize=3, markerfacecolor='None')


                        #nuv plot
                        subplots[n].errorbar(csv_f_n['t0'][good_criteria_n] - csv_f_n['t0'][0] + time_buffer, csv_f_n['mag_bgsub'][good_criteria_n],                                        yerr=[csv_f_n['mag_bgsub_err_2'][good_criteria_n], csv_f_n['mag_bgsub_err_1'][good_criteria_n]],                                         linestyle='None', marker='o', markersize=3)

                        subplots[n].errorbar(csv_f_n['t0'][bad_criteria_n] - csv_f_n['t0'][0] + time_buffer, csv_f_n['mag_bgsub'][bad_criteria_n],                                        yerr=[csv_f_n['mag_bgsub_err_2'][bad_criteria_n], csv_f_n['mag_bgsub_err_1'][bad_criteria_n]], linestyle='None',                                        marker='o', markersize=3, markerfacecolor='None')


                        #Axes
                        subplots[n].set_ylim(max_mag_limit + 0.05, min_mag_limit - 0.05)
                        subplots[n].set_xlim(0, time_buffer * 22)
                        subplots[n].set_ylabel('AB mag')
                        subplots[n].set_xlabel('Time - ' + str('%.f' % (begin_times_n[i] - time_buffer)) + ' (s)')


                        #Title
                        param = fname_n.replace('.csv', '').replace('time','').split('_')
                        subplots[n].set_title(param[0] + ' ' + param[7] + ' Visit ' + param[8] + '/' + fname_f.split('_')[7] + ' Visit ' + str(int(time_matches_n[i])) +                            ' - aper ' + param[2] + ' as, timestep ' + param[5] + (' s\nDetloc x/y %.2f/%.2f arcmin, Detrad %.2f arcmin' % (xloc, yloc, rad)) + ', Time ' + timestamp)

                    n += 1
                    if n >= numPlotsPerPage:
                        plt.tight_layout()
                        pdf.savefig(fig)
                        plt.close()
                        n = 0
                        fig, subplots = plt.subplots(numPlotsPerPage, figsize=(10,12))
                        if numPlotsPerPage == 1: subplots = [subplots] # see above comment at other plt.sbplots() call.

                # Save any plots currently left on the page.
                if n > 0:
                    # Hide any unused subplots first
                    for j in range(n, numPlotsPerPage):
                        subplots[j].axis('off')

                    # Save!
                    plt.tight_layout()
                    pdf.savefig(fig)
                    plt.close()
                    n = 0

    for csvpy in file_list_csvpy:
        try:
            plotfile(csvpy, 15)
        except KeyboardInterrupt: 
            exit()
        except:
            # Handle error
            print 'Error handling ' + csvpy + ', aperture 15: ' + str(sys.exc_info()[0]) + '(' + str(sys.exc_info()[1]) + ')'

        try:
            plotfile(csvpy, 25)
        except KeyboardInterrupt:
            exit()
        except:
            # Handle error
            print 'Error handling ' + csvpy + ', aperture 25: ' + str(sys.exc_info()[0]) + '(' + str(sys.exc_info()[1]) + ')'
                

