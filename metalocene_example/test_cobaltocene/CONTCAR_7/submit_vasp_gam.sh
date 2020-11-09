#!/bin/bash

#SBATCH -p enge
#SBATCH --qos=nogpu
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=36
#SBATCH -J /home/chr218/materials_bandgap_workflow/metalocene_example/test_cobaltocene/CONTCAR_7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=soldrivertorso@gmail.com

cd $SLURM_SUBMIT_DIR

module unload gcc/6.1.0
module load vasp/5.4.1-neb
module load mvapich2

export I_MPI_COMPATIBILITY=4
export MPICH_NO_BUFFER_ALIAS_CHECK=1
ulimit -s unlimited

srun vasp_gam > out

exit

