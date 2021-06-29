#!/bin/bash -x
#
# Script to remove data older than 2 days from GAIA's directories listed within.
#
# Author: CT Neris - 02JUN2021

dir="/home/operador/gaia/MODEL"

# Insert here the DIR NAME from where you want to delete the files
dirlist="COSMO WRF WW3_COSMO WW3_ICON WW3_GFS"

# Looping over the list removing files older than 2 days.
for dirs in $dirlist; do

	echo Removing data older than 2 days in $dir/$dirs...

	if [ $dirs == "COSMO" ]; then

		rm `find $dir/$dirs/cosmo_met5* -mtime +1 -exec ls -l {} \; | awk '{print $9}'`

	elif [ $dirs == "WRF" ]; then

		 rm `find $dir/$dirs/wrf_metarea5* -mtime +1 -exec ls -l {} \; | awk '{print $9}'`

	elif [ $dirs == "WW3_COSMO" ]; then
		 
		rm `find $dir/$dirs/ww3cosmo_met* -mtime +1 -exec ls -l {} \; | awk '{print $9}'`

	elif [ $dirs == "WW3_ICON" ]; then
		
		rm `find $dir/$dirs/ww3icon* -mtime +1 -exec ls -l {} \; | awk '{print $9}'`

	elif [ $dirs == "WW3_GFS" ]; then
		
		rm `find $dir/$dirs/ww3gfs* -mtime +1 -exec ls -l {} \; | awk '{print $9}'`

	fi
done
