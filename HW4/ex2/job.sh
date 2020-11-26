#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --exclusive
#SBATCH --partition=NODE2008

##mpirun ./a.out
mpirun --mca btl_openib_allow_ib true ./a.out < input.d
