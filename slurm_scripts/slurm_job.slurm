#!/bin/bash -l

#SBATCH --job-name=hx71    # Job name
#SBATCH --output=jobname.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=jobname.%j.err # Stderr (%j expands to jobId)
#SBATCH --ntasks=1     # Number of tasks(processes)
#SBATCH --nodes=1     # Number of nodes requested
#SBATCH --ntasks-per-node=1     # Tasks per node
#SBATCH --cpus-per-task=1     # Threads per task
#SBATCH --time=24:00:00   # walltime
#SBATCH --partition=skylake

conda run -n ADDCROSS python3 /home/agkikas/MARS_RETRIEVE_SCRIPTS/MARS-RETRIEVE-SCRIPTS/ecmwf_mars_request_pl_fc.py
