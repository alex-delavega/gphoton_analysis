# gphoton_analysis
### Alex de la Vega & Agustina Quesada (JHU)

A python module to extract and analyze time-resolved photometry from gPhoton (Million et al. 2016).

This Python 2.7 module allows users to create and analyze time-resolved GALEX UV photometry from gPhoton
using a variety of scripts and/or functions for en masse or individual object analysis. 

If you use this module in your research, please cite de la Vega, Bianchi & Quesada (in prep.). 

### BEFORE YOU BEGIN:

After cloning this repository, the following file requires user-defined environment variables:
```
gphoton_var.sh
```

**Please check this file to ensure it runs according to the parameters of your choice.**

We recommend you either `source gphoton_var.sh` wherever you clone this repo or add the following line to your `~/.bash_profile`:
```
source /path/to/gphoton_analysis/gphoton_var.sh
``` 

To make the main shell script, `create_photometry.sh`, executable, run through your Terminal (assuming you are in the gphoton_analysis directory):
```
chmod +x ./create_photometry_files.sh
```
Run it using:
```
./create_photometry_files.sh
```