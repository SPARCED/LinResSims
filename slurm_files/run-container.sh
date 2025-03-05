#!/bin/bash
###############################################################################
# Script Name:  run-container.sh
# Description:  Run the singularity container, validation 'run-all' function to
#               ensure proper function
#
# Author:       Jonah R. Huggins
# Created:      2024-12-13
# Version:      1.0
#
# Usage:        sbatch run-container.sh [options] [arguments]
#
# Options:
#
# Arguments:    -c | Number of cores for the MPI job to use during every
#                    benchmark
#                    (Example: -c 10)
#
# Requirements: slurm
#               openmpi
#
# Notes:        This file is intended for use on systems using SLURM
#               Job Scheduler.
###############################################################################
#SBATCH --job-name run-all
#SBATCH --output=run-all_output.txt # Output file
#SBATCH --error=run-all_errors.txt  # Error file
#SBATCH --nodes=2                     # Number of nodes
#SBATCH --ntasks=40                  # Total number of MPI tasks (cores)
#SBATCH --ntasks-per-node=20          # Number of MPI tasks per node, balancing load across nodes
#SBATCH --cpus-per-task=1             # Number of CPUs per task (usually 1 for MPI)
#SBATCH --mem-per-cpu=20gb            # Memory alloted to each node
#SBATCH --time 23:59:00               # Time

cd /scratch/jrhuggi/LinResSims/

module load openmpi

mpirun -n 40 singularity exec container/linressims.sif bash -c "
cd scripts
python cellpop.py --sim_config default_SPARCED.json 2>&1 | tee cellpop.log
"
