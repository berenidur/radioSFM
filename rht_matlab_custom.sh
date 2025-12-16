#!/bin/bash
#SBATCH --output="slurm/slurm-%j.out"
#SBATCH --time=86400
#SBATCH --gres=gpu:1

# Run the job
srun matlab -nosplash -nodesktop -nodisplay -r "$1; exit"