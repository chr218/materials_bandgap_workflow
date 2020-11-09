#!/bin/bash

DIR_present=$(pwd)
arr=$(find . -type d -name "CONTCAR_*")


for dir in $arr
do

	cd $dir
	~/bin/update_submissionscript_to_current_directory.sh	
	sbatch submit*
	#scancel --name $(pwd)
	cd $DIR_present
done
