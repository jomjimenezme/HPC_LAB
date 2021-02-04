#!/bin/bash
#
#SBATCH --nodes=4
#SBATCH --ntasks=25
#SBATCH --exclude=node01,node02
#SBATCH --exclusive
#SBATCH --partition=NODE2008
   mpirun  --mca btl_openib_allow_ib true ./a.out >25.d 
