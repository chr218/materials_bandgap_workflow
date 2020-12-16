#!/bin/bash

#SBATCH -p hawkmem
#SBATCH --qos=nogpu
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=52
#SBATCH -J /home/chr218/materials_bandgap_workflow/updates_11_22_2020/CONTCAR_2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=soldrivertorso@gmail.com

cd $SLURM_SUBMIT_DIR

source /etc/profile.d/zlmod.sh
module load vasp/5.4.4.pl2
ulimit -s unlimited


date > start_time

srun vasp_gam > out

date > finish_time
exit
