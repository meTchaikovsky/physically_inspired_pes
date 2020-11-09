#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=debug2
#SBATCH --nodes=1

Host=`scontrol show job $SLURM_JOB_ID | grep ' NodeList' | awk -F'=' '{print $2}'`
python create_db.py $Host
