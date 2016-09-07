# gphoton_analysis
### Alex de la Vega & Halley Cromley (JHU)

A python module to analyze output from gPhoton (Million et al. 2016).

This Python 2.7 module allows users to create and analyze output files from gPhoton
using a variety of scripts and/or functions for en masse or individual object analysis. 

### BEFORE YOU BEGIN:

The following scripts require user-defined paths to save output files:
```
gfind_to_times.py
times_to_csv.py
create_photometry_files.sh
stat_plot_csv.py
```

**Please check these files to ensure they run and save to the paths you desire.**

To make the .sh file executable, run through your Terminal (assuming you are in the gphoton_analysis directory):
```
chmod +x ./create_photometry_files.sh
```
Run it using:
```
./create_photometry_files.sh
```

**Make sure to log progress in log directory!**
