#!/bin/bash
#
#SBATCH --nodes=8
#SBATCH --ntasks=64
#SBATCH --exclusive
#SBATCH --partition=NODE2008

##mpirun ./a.out
mpirun --mca btl_openib_allow_ib true ./a.out < input.d
