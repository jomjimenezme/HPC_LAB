#!/bin/bash
#
#SBATCH --nodes=8
#SBATCH --exclude=node05
#SBATCH --exclusive
#SBATCH --partition=NODE2008
i=4
while [ $i -lt 65 ]
do
   mpirun -np $i --mca btl_openib_allow_ib true ./scalingjacobi  >scaling$i.d
   i=$[$i*4]
   wait
done
