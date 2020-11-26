#!/bin/bash
#
#SBATCH --nodes=2
#SBATCH --account=training2008
for i in {2..64}
do
   mpirun -np $i ./getdata  > getdata$i.out
done
