"""
This file reads in a CSV of targets and prepares a python script, 
which calls gFind and saves the output to a directory where
time ranges for each target are stored.

Alexander de la Vega -- 5 / 9 / 2016
"""

import os # to change directories
import csv # to read in files
from astropy.io import ascii

time_path = './times/' # path to store time ranges for each target
if not os.path.exists(time_path): # check whether times directory exists
    os.makedirs(time_path) # if not, create directory

find_path = './find/' # path to store gFind output for each target
if not os.path.exists(find_path): # check whether gFind output directory exists
    os.makedirs(find_path) # if not, create directory

csv_path = '' # path to CSV with targets
csv_f = ascii.read(csv_path)

for i in range(len(csv_f)): # ignoring the first line
    objid = str(csv_f['objid'][i]) # name for files
    ra, dec = str(csv_f['ra'][i]), str(csv_f['dec'][i]) # right ascension and declination
    fuvmag, nuvmag = str(csv_f['fuv_mag'][i]), str(csv_f['nuv_mag'][i])
    
    g = open(find_path + objid + '_find.py', 'w') # open file for gFind 
    g.write('# RA: ' + ra + '\n') 
    g.write('# Dec: ' + dec + '\n')
    g.write('# FUV mag: ' + fuvmag + '\n')
    g.write('# NUV mag: ' + nuvmag + '\n\n')
    g.write('import gPhoton\n') # importing gPhoton
    g.write('import sys \n') # capture output from shell
    g.write('from cStringIO import StringIO \n\n')
    g.write('# setup the environment \n'\
            'backup = sys.stdout \n\n')
    g.write('# #### \n'\
            'sys.stdout = StringIO()     # capture output \n\n')

    g.write("gPhoton.gFind(band='FUV',skypos=[" + ra + "," + dec + "])\n\n") # find target observation times - FUV
    
    g.write("gPhoton.gFind(band='NUV',skypos=[" + ra + "," + dec + "])\n\n") # " " - NUV
    
    # below records output from shell in a new python file and sends it to /times directory

    g.write('out = sys.stdout.getvalue() # release output \n\n'\

            'sys.stdout.close()  # close the stream \n'\
            'sys.stdout = backup # restore original stdout \n\n'\

            'print out.upper()   # post processing \n \n'\

            "h = open('" + time_path + objid + "_times.py', 'w') \n\n"\

            'h.write(out.upper()) \n\n'\

            'h.close()')
    
    g.close()