#!/bin/bash
#SBATCH --job-name=grabflucts_array
#SBATCH --partition=sixhour
#SBATCH --output=grabflucts_array%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --constraint=intel
#SBATCH --mem=30G
#SBATCH --time=06:00:00
#SBATCH --array=0-6

module load Dir_Calc_Fluct 

cd $SLURM_SUBMIT_DIR

echo Time is `date` > array_$SLURM_ARRAY_TASK_ID.o
echo Directory is `pwd` >> array_$SLURM_ARRAY_TASK_ID.o

grab_flucts.py flucts.inp $SLURM_ARRAY_TASK_ID
