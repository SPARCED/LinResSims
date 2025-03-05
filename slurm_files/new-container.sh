#!/bin/bash
###############################################################################
# Script Name:  new-container.sh
# Description:  Build the singularity definition file into a container for use.
#
# Author:       Jonah R. Huggins
# Created:      2024-08-29
# Version:      1.0
#
# Usage:        sbatch new-container.sh [options] [arguments]
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
#SBATCH --job-name container-build
#SBATCH --output=build_output.txt # Output file
#SBATCH --error=build_errors.txt  # Error file
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Total number of MPI tasks (cores)
#SBATCH --ntasks-per-node=1           # Number of MPI tasks per node, balancing load across nodes
#SBATCH --cpus-per-task=1             # Number of CPUs per task (usually 1 for MPI)
#SBATCH --mem-per-cpu=20gb            # Memory alloted to each node
#SBATCH --time 04:59:00               # Time

cd /scratch/jrhuggi/LinResSims/
singularity build --fakeroot container/linressims.sif container/linressim.def
