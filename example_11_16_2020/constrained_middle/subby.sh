#!/bin/bash

DIR_present=$(pwd)
arr=$(find . -type d -name "CONTCAR_*")


for dir in $arr
do

	cd $dir
	#sed -i "s/engc/enge/g" /submit_vasp_gam.sh
	cp $DIR_present/submit_vasp_gam.sh .
	~/bin/update_submissionscript_to_current_directory.sh	
	sbatch /submit_vasp_gam.sh
	#scancel --name $(pwd)
	cd $DIR_present
done
