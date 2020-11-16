#!/bin/bash

DIR_present=$(pwd)
arr=$(find . -type d -name "CONTCAR_*")


for dir in $arr
do

	cd $dir
	cp $DIR_present/submit_vasp_gam.sh .	
	~/bin/update_submissionscript_to_current_directory.sh	
	sbatch submit_vasp_gam.sh
	#scancel --name $(pwd)
	cd $DIR_present
done
