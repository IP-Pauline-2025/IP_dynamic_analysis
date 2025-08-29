#!/bin/bash
#SBATCH --job-name=wind-dyn
#SBATCH --mail-user=bittner@irz.uni-hannover.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=uq_model_%j.out     # Output file
#SBATCH --error=uq_model_%j.err      # Error file
#SBATCH --ntasks=40                  # Number of tasks (adjust as needed)

# Navigate to the working directory
cd $SLURM_SUBMIT_DIR

module load OpenSees

# Run Julia script with parallel execution
julia --project=. --threads=$SLURM_NTASKS dynamic_analysis_2.jl