#!/bin/bash
#
#SBATCH --nodes=8
#SBATCH --exclude=node05
#SBATCH --exclusive
#SBATCH --partition=NODE2008
i=2
while [ $i -lt 65 ]
do
   mpirun -np $i --mca btl_openib_allow_ib true ./a.out  >timing$i.out
   i=$[$i*2]
   wait
done
