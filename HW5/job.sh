#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --exclude=node05
#SBATCH --exclusive
#SBATCH --partition=NODE2008
   mpirun -np 1 --mca btl_openib_allow_ib true ./a.out  >timing1.out
