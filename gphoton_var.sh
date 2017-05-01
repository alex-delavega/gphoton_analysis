#!/bin/bash

# GPHOTON_VAR.SH
# This file lists the main parameters necessary for extracting time-resolved photometry
# from gPhoton. This includes: main directories, aperture and background annulus radii, 
# and time bin size. 
# Alexander de la Vega & Agustina Quesada
# last modified 1 May 2017

# below are main directories 
export GPHOTON_HOME="/Users/Alex/Desktop/new_gphoton/gphoton_analysis/"
export GPHOTON_INPUT="/Users/Alex/Desktop/new_gphoton/gphoton_analysis/input/"
export GPHOTON_OUTPUT="/Users/Alex/Desktop/new_gphoton/gphoton_analysis/output/"
export GPHOTON_FIND="/Users/Alex/Desktop/new_gphoton/gphoton_analysis/find/"
export GPHOTON_TIMES="/Users/Alex/Desktop/new_gphoton/gphoton_analysis/times/"
export GPHOTON_CSVPY="/Users/Alex/Desktop/new_gphoton/gphoton_analysis/CSVpy/"
export GPHOTON_PLOTS="/Users/Alex/Desktop/new_gphoton/gphoton_analysis/plots/"

# below are parameters for gPhoton photometry
# time bin size IN SECONDS
export GPHOTON_TIME_BIN=5

# aperture radius IN ARCSECONDS
export GPHOTON_APERTURE=15

# background inner radius IN ARCSECONDS
export GPHOTON_BCKG_INNER=30

# background outer radius IN ARCSECONDS
export GPHOTON_BCKG_OUTER=45

# below are specific parameters for gPhoton input files
# the input file of choice
export GPHOTON_INPUT_FILE="${GPHOTON_INPUT}/targets.csv"

# object id keyword in input file
export GPHOTON_OBJID="objid"

# right ascension keyword in input file
export GPHOTON_RA="ra"

# declination keyword in input file
export GPHOTON_DEC="dec"