#!/bin/bash
#
#SBATCH --nodes=2
#SBATCH --account=training2008
for i in {2..64}
do
   mpirun -np $i ./whole  >whole$i.out
done
