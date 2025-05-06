#!/bin/bash -l

#SBATCH -J S2_fit_models
#SBATCH -A naiss2024-22-1028  # Replace with your PDC project allocation
#SBATCH -t 24:00:00          # Adjust time as needed
#SBATCH --nodes=1            # Start with 1 node, increase if needed
#SBATCH --ntasks-per-node=1  # For Dardel thin nodes (adjust based on node type)
#SBATCH -o output_%a.txt
#SBATCH -e error_%a.txt
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH  -p main
#SBATCH --cpus-per-task=5
#SBATCH --array=1-3 # remove for running models one by one

# --- Load necessary modules (if any) ---
#  module load ...  # Add any modules needed (R, etc.)
module load PDC/23.12 R/4.4.1-cpeGNU-23.12

# --- Define Variables ---
output_dir="/cfs/klemming/projects/supr/naiss2024-22-1028/models/" 
unfitted_models_file="/cfs/klemming/projects/supr/naiss2024-22-1028/models/unfitted_models.RData"
nParallel=4
nChains=4

# --- Calculate model index based on Slurm array task ID ---
model_index=$SLURM_ARRAY_TASK_ID # remove for running models one by one

# --- Run the R script ---
Rscript "/cfs/klemming/projects/supr/naiss2024-22-1028/S2_fit_models_hpc.R" \
    "$output_dir" \
    "$unfitted_models_file" \
    "$nChains" \
    "$nParallel" \
    "$model_index" # remove for running models one by one

  