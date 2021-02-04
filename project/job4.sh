#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --exclude=node01,node02
#SBATCH --exclusive
#SBATCH --partition=NODE2008
   mpirun  --mca btl_openib_allow_ib true ./a.out >4.d 
