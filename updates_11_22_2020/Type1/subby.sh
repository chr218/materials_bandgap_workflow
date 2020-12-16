#!/bin/bash

DIR_present=$(pwd)
arr=$(find . -type d)


for dir in $arr
do

	cd $dir
	cp ~/bin/submit_vasp_gam_HAWK.sh .
	~/bin/update_submissionscript_to_current_directory.sh
	cp $DIR_present/INCAR .
	sbatch submit_vasp_gam_HAWK.sh
	cd $DIR_present

done
