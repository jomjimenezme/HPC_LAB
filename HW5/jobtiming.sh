#!/bin/bash
#
#SBATCH --nodes=8
#SBATCH --exclusive
#SBATCH --partition=NODE2008
for i in {2..64}
do
   mpirun -np $i --mca btl_openib_allow_ib true ./a.out  >timing$i.out
done
