#!/bin/bash
#SBATCH --job-name=hex_k-6
#SBATCH --partition=reserved
#SBATCH --account=rac-2018-hpcg1742
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pavithran.sridhar@gmail.com
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8192
#SBATCH --array=0-5
echo "Running task" $SLURM_ARRAY_TASK_ID
./instances.sh $SLURM_ARRAY_TASK_ID