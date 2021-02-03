#!/bin/bash
#
#SBATCH --nodes=8
#SBATCH --exclude=node05
#SBATCH --exclusive
#SBATCH --partition=NODE2008
i=2
while [ $i -lt 9 ]
do
  mpirun -np $[$i*$i] --mca btl_openib_allow_ib true ./strongscaling_SOR  >srongscaling_SOR$[$i*$i].d
  i=$[$i+1]
  wait
done

