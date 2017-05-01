
#!/bin/bash
# this bash file executes python files to: 
# 1) extract info from target list CSV to find observation times for targets via gFind
# 2) saves gFind output in times folder
# 3) uses observation time data to produce CSVs through gAperture
# 4) final ready-to-analyze CSVs in CSV folder

# Alexander de la Vega & Agustina Quesada
# last modified: 1 May 2017

cd $GPHOTON_HOME

python csv_to_gfind.py # execute initial file that launches gFind

cd $GPHOTON_FIND # go to /find folder 

for f in $(ls); do # activate each python script to send gFind output to /times folder
	python $f; 
done

cd $GPHOTON_HOME # go back to main directory

python times_to_csv.py # execute python script that prepares files that will produce CSVs

cd $GPHOTON_CSVPY # go to directory with scripts that produce CSVs

for f in $(ls); do # run script to make CSVs
	python $f; 
done