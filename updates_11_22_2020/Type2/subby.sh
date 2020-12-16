#!/bin/bash

DIR_present=$(pwd)
arr=$(find . -type d)


for dir in $arr
do

	cd $dir
	~/bin/update_submissionscript_to_current_directory.sh
	sbatch submit_vasp_gam.sh
	cd $DIR_present

done
