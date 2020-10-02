#!/bin/bash
#SBATCH --job-name=water_jump
#SBATCH --partition=thompson,laird,shontz
#SBATCH --constraint=intel
#SBATCH --output=output.log
#SBATCH --mail-user=piskuliche@ku.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=50:00:00


module purge
module load compiler/intel/18 intel-mpi/18
module load lammps/12Dec2018
echo Time is `date`
echo Directory is `pwd`

echo "Running on $SLURM_JOB_NODELIST nodes using $SLURM_CPUS_ON_NODE cores on each node"

mpirun lmp_mpi -in in.water > out

cd ../
echo Ending Time is `date`
exit 0
