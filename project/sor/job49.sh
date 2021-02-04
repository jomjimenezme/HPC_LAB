#!/bin/bash
#
#SBATCH --nodes=7
#SBATCH --ntasks=49
#SBATCH --exclude=node01,node02
#SBATCH --exclusive
#SBATCH --partition=NODE2008
   mpirun  --mca btl_openib_allow_ib true ./a.out >49.d 
