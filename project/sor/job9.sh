#!/bin/bash
#
#SBATCH --nodes=2
#SBATCH --ntasks=9
#SBATCH --exclude=node01,node02
#SBATCH --exclusive
#SBATCH --partition=NODE2008
   mpirun  --mca btl_openib_allow_ib true ./a.out > 9.d 
