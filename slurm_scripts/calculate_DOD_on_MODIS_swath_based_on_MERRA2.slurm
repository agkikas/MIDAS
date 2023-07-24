#!/bin/bash -l

#SBATCH --job-name=1A2022    # Job name
#SBATCH --output=jobname.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=jobname.%j.err # Stderr (%j expands to jobId)
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=agkikas@noa.gr
#SBATCH --ntasks=1     # Number of tasks(processes)
#SBATCH --nodes=1     # Number of nodes requested
#SBATCH --ntasks-per-node=1     # Tasks per node
#SBATCH --cpus-per-task=1     # Threads per task
#SBATCH --mem=6G #total memory per node
#SBATCH --time=1-23:00:00   # walltime
#SBATCH --partition=cascadelake

conda run -n MIDAS python3 /home/agkikas/optimized_MIDAS_scripts/calculate_DOD_on_MODIS_swath_based_on_MERRA2.py
