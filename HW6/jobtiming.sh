#!/bin/bash
#
#SBATCH --nodes=8
#SBATCH --exclude=node05
#SBATCH --exclusive
#SBATCH --partition=NODE2008
mpirun -np 1 --mca btl_openib_allow_ib true ./timing  >timing1.out
i=2
while [ $i -lt 65 ]
do
   mpirun -np $i --mca btl_openib_allow_ib true ./timing  >timing$i.out
   i=$[$i*2]
   wait
done
