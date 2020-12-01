#!/bin/bash
#
#SBATCH --nodes=8
#SBATCH --exclusive
#SBATCH --exclude=node05
#SBATCH --partition=NODE2008

##mpirun ./a.out
for i in {2..64}
do
   mpirun --mca btl_openib_allow_ib true -np $i ./reduce  >reduce$i.out
done
