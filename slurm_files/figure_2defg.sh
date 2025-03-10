#!/bin/bash
###############################################################################
# Script Name:  figure_2defg.sh
# Description:  Runs all necessary simulations for reproducing drug-dose response
#               curves for Figure 2D-G. 
#
# Author:       Jonah R. Huggins
# Created:      2025-03-06
# Version:      1.0
#
# Usage:        sbatch figure_2defg.sh [options] [arguments]
#
# Options:
#
# Arguments:    
#
# Requirements: slurm
#               openmpi
#
# Notes:        This file is intended for use on systems using SLURM
#               Job Scheduler.
###############################################################################
#SBATCH --job-name fig_2defg
#SBATCH --output=sim_output.txt # Output file
#SBATCH --error=sim_errors.txt  # Error file
#SBATCH --nodes=2                     # Number of nodes
#SBATCH --ntasks=40                    # Total number of MPI tasks (cores)
#SBATCH --ntasks-per-node=1           # Number of MPI tasks per node, balancing load across nodes
#SBATCH --cpus-per-task=1             # Number of CPUs per task (usually 1 for MPI)
#SBATCH --mem-per-cpu=20gb            # Memory alloted to each node
#SBATCH --time 71:59:00               # Time

# load necessary modules
module load openmpi

# Move to the appropriate directory (User specific)
cd /scratch/jrhuggi/LinResSims

#---------------------Drugs & Concentrations Variables-------------------------#
drugs=("alpel_EC" "nerat_EC" "trame_EC" "palbo_EC")

concentrations=(0.0 0.001 0.003162 0.01 0.031623
                0.1 0.316228 1.0 3.162278 10.0)

#------------------------Main Operation----------------------------------------#
for drug in "${drugs[@]}"; 
do
    for conc in "${concentrations[@]}"; 
    do       
        # set the current drug and dose concentration.
        python scripts/update_json.py -p sim_configs/drs_SPARCED.json --drug "$drug" --dose "$conc"

        mpirun -n 100 singularity exec container/linressims.sif bash -c "
        cd scripts
        python cellpop.py --sim_config drs_SPARCED.json 
        "
    done
done