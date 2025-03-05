#!/bin/bash
###############################################################################
# Script Name:  compile-container.sh
# Description:  Compile the SPARCED model for use.
#
# Author:       Jonah R. Huggins
# Created:      2024-08-29
# Version:      1.0
#
# Usage:        sbatch compile-container.sh [options] [arguments]
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
#SBATCH --job-name compile-model
#SBATCH --output=compile_output.txt # Output file
#SBATCH --error=compile_errors.txt  # Error file
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Total number of MPI tasks (cores)
#SBATCH --ntasks-per-node=1           # Number of MPI tasks per node, balancing load across nodes
#SBATCH --cpus-per-task=1             # Number of CPUs per task (usually 1 for MPI)
#SBATCH --mem-per-cpu=20gb            # Memory alloted to each node
#SBATCH --time 04:59:00               # Time

cd /scratch/jrhuggi/LinResSims/
singularity exec container/linressims.sif bash -c "
cd scripts
python createModel.py 2>&1 | tee createModel.log
"
