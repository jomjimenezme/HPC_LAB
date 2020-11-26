#!/bin/bash
#
#SBATCH --nodes=3
#SBATCH --ntasks=24
#SBATCH --exclusive
#SBATCH --partition=NODE2008

##mpirun ./a.out
mpirun --mca btl_openib_allow_ib true ./a.out 
