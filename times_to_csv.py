"""
This script reads recorded output from gFind (observation times) and formats input to create 
another python file, which launches gAperture to produce CSVs for each observation time.

Alexander de la Vega -- 5 / 9 / 2016
"""

import os # to change directories
import csv # to read in files
from astropy.io import ascii

time_path = './times/' # path to store time ranges for each target
if not os.path.exists(time_path): # check whether times directory exists
    os.makedirs(time_path) # if not, create directory

files = os.listdir(time_path) # list of files containing observation time output from gFind
file_len = len(files) # number of files

find_path = './find/' # path to store gFind output for each target
if not os.path.exists(find_path): # check whether gFind output directory exists
    os.makedirs(find_path) # if not, create directory

csvpy_path = './CSVpy/' 
if not os.path.exists(csvpy_path): # check whether directory for scripts to create photometry CSVs exists
    os.makedirs(csvpy_path) # if not, create directory

out_path = './out/'
if not os.path.exists(out_path): # check whether output directory exists
    os.makedirs(out_path) # if not, create directory

aper = 0.00416667 # aperture radius in deg. 
stepsize = 5.0 # time step in sec

apstring = '15' # aperture radius in arcsec
stepstr = '5' # time step in sec

inner = 0.00833333 # inner radius of annulus in deg.
outer = 0.0125 # outer radius of annulus in deg. 

for i in range(file_len): # go through each file
    if files[i].endswith('times.py'): # times.py signifies gFind output
        times = open(time_path + files[i]) # observation time file
        content = times.readlines() # read in the file
        lines = len(content) # number of lines in file
        
        for l in range(len(content)): # for each line in file
            if content[l][0] == 'F': # if the very first character is 'F', then there are FUV observation times
                f_index = l # index where the first letter is 'F'
                f_numexp = int(content[l].split()[4]) # .split() separates the string into non-white space characters; number of exposures is the fifth element in the split string
    
            if content[l][0] == 'N': # same as above (SAA) but for NUV
                n_index = l
                n_numexp = int(content[l].split()[4])
        
        name = files[i][:-9] # essentially the objid from csv_to_gfind.py - removes the 'times.py' at the end
        
        radecf = open(find_path + name + '_find.py') # open find.py file to extract RA & Dec
        coord_info = radecf.readlines() # read file
        
        ra = coord_info[0].split()[2] # the right ascension
        dec = coord_info[1].split()[2] # SAA for declination
        #fuvmag = coord_info[2].split()[3] # '' photometric bands
        #nuvmag = coord_info[2].split()[3]
        
        radecf.close() # close file
    
        csvf = open(csvpy_path + name + '_csv.py', 'w') # file that gets launched to activate gAperture
        
        csvf.write('# Object id: ' + name +'\n') # header
        csvf.write('# RA: '+ ra +'\n# Dec: '+ dec +'\n\n') # header
        csvf.write('#FUV exposure times: ' + str(f_numexp) +'\n') # list of FUV exposure times
    
        for a in range(f_index + 1, f_numexp + 1): # print out time ranges
            csvf.write('FUV_t'+ str(a) +'= ['+ content[a].split()[1].replace(',', '') + ', ' + content[a].split()[2].replace(',', '') + '] \n')
        
        csvf.write('\n')
        csvf.write('#NUV exposure times: '+ str(n_numexp) +'\n') # same as above but for NUV
        
        for b in range(n_index + 1, lines):
            csvf.write('NUV_t'+ str(b - n_index) +'= ['+ content[b].split()[1].replace(',', '') + ', ' + content[b].split()[2].replace(',', '') + '] \n')
        
        csvf.write('\n\n')
        csvf.write('import gPhoton \n\n') # import gAperture
        csvf.write('def main(): \n') # name of dummy function to be called at end of file
        
        f_num_exists = 'f_numexp' in locals() or 'f_numexp' in globals() # if target has FUV data
        n_num_exists = 'n_numexp' in locals() or 'n_numexp' in globals() # if target has NUV data

        if f_num_exists: # provided FUV data
            for c in range(f_numexp): # print a gAperture call for each FUV observation time
                csvf.write("\tgPhoton.gAperture(band='FUV', "\
                        "skypos=[" + str(ra) + "," + str(dec) +"], "\
                        "stepsz=" + str(stepsize) + ", "\
                        "csvfile=" + out_path + name + "_aper_" + apstring + "_as_timestep_" + stepstr + "_sec_FUV_time" + str(c+1) + ".csv', "\
                        "radius=" + str(aper) + ", annulus=[" + str(inner) + "," + str(outer) + "], trange=FUV_t" + str(c+1) +")\n")
       
        csvf.write('\n') # newline
        
        if n_num_exists: # same as above but for NUV
            for d in range(n_numexp):
                csvf.write("\tgPhoton.gAperture(band='NUV', "\
                        "skypos=[" + str(ra) + "," + str(dec) + "], "\
                        "stepsz=" + str(stepsize) + ", "\
                        "csvfile=" + out_path + name + "_aper_" + apstring + "_as_timestep_" + stepstr + "_sec_NUV_time" + str(d+1) + ".csv', "\
                        "radius=" + str(aper) + ", annulus=[" + str(inner) + "," + str(outer) + "], trange=NUV_t" + str(d+1) + ")\n")

        csvf.write('\n')
        csvf.write("if __name__ == '__main__':\n") # execute file
        csvf.write('\tmain()')
        
        times.close() # close observation time file
	csvf.close() # close new python file