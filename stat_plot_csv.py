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

Alexander de la Vega (JHU)
Latest version: 2 November 2016
"""

from __future__ import division # python 2.x usually has issues dividing
import numpy as np # for computation / array formatting
import os # move through directories
import sys # read in command line arguments
#from operator import itemgetter # used to retrieve elements from arrays (look in plots section)

import useful_funcs as funcs # functions used throughout gPhoton analysis
from astropy.io import ascii

def is_float(x): # check whether values from CSVs are floats or not
    try:
        float(x)
        return True
    except ValueError:
        return False

time_path = './times/' # directory with files that contain observation times
if not os.path.exists(time_path): # check whether times directory exists
    print 'Please check input for times path'
    sys.exit(0)

files = os.listdir(time_path) # list of all files in times directory (i.e. each object of interest)
num_file = len(files) # number of files in times directory

csv_path = '../new_csv/' # path that holds photometry CSVs
if not os.path.exists(csv_path): # check whether CSV directory exists
    print 'Please check input for CSV path'
    sys.exit(0)
file_list = os.listdir(csv_path)

files_csv = os.listdir(csv_path) # list of CSVs in CSV directory
num_csv = len(files_csv) # number of CSVs

out_path = './output/' # output path
if not os.path.exists(out_path): # check whether output directory exists
    os.makedirs(out_path) # if not, create directory

find_path = './find/' # path of find.py files - includes ra, dec, mag
if not os.path.exists(find_path): # check whether gFind output directory exists
    print 'Please check input for find path'
    sys.exit(0)

if len(sys.argv) > 2: # if input is given as .txt file
    if sys.argv[2].endswith('.txt'):
        inputf = open(sys.argv[2], 'r')
        input_cont = inputf.readlines()
        num_file = len(input_cont)
        file_list = [n.replace('\n', '') for n in input_cont]
        inputf.close()

aperture = 15.0 # aperture radius in arcsec
time_step = 5.0 # time step in sec

aper = '15' # vector of aperture radii strings in arcsec; used to construct file names for output files
ts = '5' # time step in sec - vector of strings for file name construction mentioned above
band = ['FUV', 'NUV'] # vector of photometric bands - for aforementioned file name construction

# if user passes 'stats' to the command line - produce stats output .txt file  
if sys.argv[1] == 'stats':
    from scipy.signal import argrelextrema
    from astropy.convolution import convolve_fft, Box1DKernel

    kernel = Box1DKernel(3)
    char_freq = 0.00882353

    good_flag = '0.0'

    stats_file = open(out_path + 'fft_stats.csv', 'w')
    stats_file.write('objid,band,visit,duration,min_mag,max_mag,mean_mag_bgsub,max_mag_var,max_mag_var_sigma,flag0_max_mag_var,flag0_max_mag_var_sigma,max_counts_var_ratio,max_counts_var_err_ratio,' +\
                 'mean_det_rad,err_det_rad,num_harmonics,num_other_peaks,flags,frac_bins_flags,tot_frac_bins_flags\n')

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
            max_counts_var_err_ratio = csv_f['cps_err'][good_exp][max_counts_index] *\
            csv_f['exptime'][good_exp][max_counts_index] / mean_counts
            
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
    from pylab import rcParams
    rcParams['figure.figsize'] = 8, 6

    min_mag_limit = 20
    max_mag_limit = 5

    csvpy_path = './CSVpy/'
    if not os.path.exists(csvpy_path):
        print 'Please check input for CSVpy path'
        sys.exit(0)

    file_list_csvpy = ['2412102149450243133_csv.py'] #os.listdir(csvpy_path)

    for f in file_list_csvpy:
        if f.endswith('csv.py'):
            csvpy_f = f
            csvpy_file = open(csvpy_path + f, 'r')
            content = csvpy_file.readlines()
            csvpy_file.close()

            for i in range(len(content)):
                if len(content[i]) > 1:
                    if content[i].split()[0] == '#FUV':
                        f_line = i
                    if content[i].split()[0] == '#NUV':
                        n_line = i
                        break

            f_exp, n_exp = int(content[f_line].split()[3]), int(content[n_line].split()[3])

            begin_times_f, end_times_f = np.empty(f_exp), np.empty(f_exp)
            begin_times_n, end_times_n = np.empty(n_exp), np.empty(n_exp)

            for i in range(f_line + 1, f_line + f_exp + 1):
                begin_times_f[i - f_line - 1] = float(content[i].replace('[','').replace(']','').replace(',','').split()[1])
                end_times_f[i - f_line - 1] = float(content[i].replace('[','').replace(']','').split()[2])
    
            for j in range(n_line + 1, n_line + n_exp + 1):
                begin_times_n[j - n_line - 1] = float(content[j].replace('[','').replace(']','').replace(',','').split()[1])
                end_times_n[j - n_line - 1] = float(content[j].replace('[','').replace(']','').split()[2])

            time_matches_f = np.asarray([i for i in range(1, f_exp + 1)])
    
            if n_exp >= f_exp:
                time_matches_n = np.zeros(n_exp)
                count = 0
                for k in range(n_exp):
                    if begin_times_n[k] in begin_times_f:
                        count += 1
                        time_matches_n[k] = count
            else:
                print 'Need to investigate'

            for f in file_list:
                if f.startswith(csvpy_f.split('_')[0]):
                    csv_f = ascii.read(csv_path + f)
                    f_min_mag = np.min(csv_f['mag_bgsub'] - csv_f['mag_bgsub_err_1'])
                    f_max_mag = np.max(csv_f['mag_bgsub'] + csv_f['mag_bgsub_err_2'])
                    if f_min_mag < min_mag_limit:
                        min_mag_limit = f_min_mag
                    if f_max_mag > max_mag_limit:
                        max_mag_limit = f_max_mag

            if n_exp >= f_exp:
                for i in range(n_exp):
                    if time_matches_n[i] == 0:
                        fname = csvpy_f.split('_')[0] + '_aper_15_as_timestep_5_sec_NUV_time' + str(i + 1) + '.csv'
                        #file_name = './new_csv/' + fname
                        csv_f = ascii.read(csv_path + fname)
            
                        time_range = csv_f['t0'][len(csv_f) - 1] - csv_f['t0'][0]
                        time_buffer = 0.05 * time_range
                        flags = funcs.find_unique_elements(csv_f['flags'])
            
                        good_criteria = (csv_f['exptime'] >= .75 * 5)
                        bad_criteria = (csv_f['exptime'] < .75 * 5)
            
                        if time_range <= 200: # if the length of the visit <= 200 sec
                            minorLocator = MultipleLocator(5) # set minor tick marks every 5 sec
                        else:
                            minorLocator = MultipleLocator(20) # set minor tick marks every 20 sec
                
                        fig, ax = plt.subplots(figsize=(6,4))
                        ax.xaxis.set_minor_locator(minorLocator) # plot with minor tick marks

                        param = fname.replace('.csv', '').replace('time','').split('_')
            
                        colors = plt.cm.rainbow(np.linspace(0, 1, len(flags)))
            
                        for j in range(len(flags)):
                            plt.errorbar(csv_f['t0'][csv_f['flags'] == flags[j]] - csv_f['t0'][0] + time_buffer, csv_f['mag_bgsub'][csv_f['flags'] == flags[j]],\
                                        yerr=[csv_f['mag_bgsub_err_2'][csv_f['flags'] == flags[j]], csv_f['mag_bgsub_err_1'][csv_f['flags'] == flags[j]]], linestyle='None',\
                                        marker='o', markersize=3, color=colors[j], label=str(flags[j]).replace('.0', ''))
            
                        plt.ylim(max_mag_limit + 0.05, min_mag_limit - 0.05)
                        plt.xlim(0, time_buffer * 22)
                        plt.ylabel('AB mag')
                        plt.xlabel('Time - ' + str('%.f' % (csv_f['t0'][0] - time_buffer)) + ' (s)')
                        plt.legend(loc='best', numpoints=1, title='Flags')
                        plt.title(param[0] + ' ' + param[7] + ' Visit ' + param[8] +\
                            ' - aper ' + param[2] + ' as, timestep ' + param[5] + ' s')
                        fig_name = out_path + param[0] + '_ts_' + ts + '_aper_' + aper + '_as_visit' + param[8] + '.pdf' # name of plot file
                        plt.savefig(fig_name, format='pdf', bbox_inches='tight', dpi=150) # save plot
                        plt.clf()
                        plt.close('all')
                    else:
                        fname_n = csvpy_f.split('_')[0] + '_aper_15_as_timestep_5_sec_NUV_time' +\
                            str(i + 1) + '.csv'
                        fname_f = csvpy_f.split('_')[0] + '_aper_15_as_timestep_5_sec_FUV_time' +\
                            str(int(time_matches_n[i])) + '.csv'
                        csv_f_n = ascii.read(csv_path + fname_n)
                        csv_f_f = ascii.read(csv_path + fname_f)
            
                        #f_index = int(time_matches_n[i]) - 1
            
                        if len(csv_f_f) >= len(csv_f_n):
                            time_range = csv_f_f['t0'][len(csv_f_f) - 1] - csv_f_f['t0'][0]
                        else:
                            time_range = csv_f_n['t0'][len(csv_f_n) - 1] - csv_f_n['t0'][0]
                        time_buffer = 0.05 * time_range
            
                        good_criteria = (csv_f['exptime'] >= .75 * 5)
                        bad_criteria = (csv_f['exptime'] < .75 * 5)
            
                        if time_range <= 200: # if the length of the visit <= 200 sec
                            minorLocator = MultipleLocator(5) # set minor tick marks every 5 sec
                        else:
                            minorLocator = MultipleLocator(20) # set minor tick marks every 20 sec
                
                        fig, ax = plt.subplots(figsize=(6,4))
                        ax.xaxis.set_minor_locator(minorLocator) # plot with minor tick marks

                        param = fname_n.replace('.csv', '').replace('time','').split('_')
            
                        f_flags = funcs.find_unique_elements(csv_f_f['flags'])
                        f_colors = plt.cm.rainbow(np.linspace(0, 1, len(f_flags)))
            
                        for j in range(len(f_flags)):
                            plt.errorbar(csv_f_f['t0'][csv_f_f['flags'] == f_flags[j]] - csv_f_f['t0'][0] + time_buffer, csv_f_f['mag_bgsub'][csv_f_f['flags'] == f_flags[j]],\
                                        yerr=[csv_f_f['mag_bgsub_err_1'][csv_f_f['flags'] == f_flags[j]], csv_f_f['mag_bgsub_err_1'][csv_f_f['flags'] == f_flags[j]]], linestyle='None',\
                                        marker='o', markersize=3, color=f_colors[j], label=str(f_flags[j]).replace('.0', ''))
            
                        n_flags = funcs.find_unique_elements(csv_f_n['flags'])
                        n_colors = plt.cm.rainbow(np.linspace(0, 1, len(n_flags)))    
            
                        for j in range(len(n_flags)):
                            plt.errorbar(csv_f_n['t0'][csv_f_n['flags'] == n_flags[j]] - csv_f_n['t0'][0] + time_buffer, csv_f_n['mag_bgsub'][csv_f_n['flags'] == n_flags[j]],\
                                        yerr=[csv_f_n['mag_bgsub_err_1'][csv_f_n['flags'] == n_flags[j]], csv_f_n['mag_bgsub_err_1'][csv_f_n['flags'] == n_flags[j]]], linestyle='None',\
                                        marker='^', markersize=3, color=n_colors[j], label=str(n_flags[j]).replace('.0', ''))
            
                        plt.ylim(max_mag_limit + 0.05, min_mag_limit - 0.05)
                        plt.xlim(0, time_buffer * 22)
                        plt.ylabel('AB mag')
                        plt.xlabel('Time - ' + str('%.f' % (begin_times_n[i] - time_buffer)) + ' (s)')
                        plt.legend(loc='best', numpoints=1, title='Flags')
                        plt.title(param[0] + ' ' + param[7] + ' Visit ' + param[8] + '/' + fname_f.split('_')[7] + ' Visit ' + str(int(time_matches_n[i])) +\
                            '- aper ' + param[2] + ' as, timestep ' + param[5] + ' s')
                        fig_name = out_path + param[0] + '_ts_' + ts + '_aper_' + aper + '_as_visit' + param[8] + '.pdf' # name of plot file
                        plt.savefig(fig_name, format='pdf', bbox_inches='tight', dpi=150) # save plot
                        plt.clf()
                        plt.close('all')



