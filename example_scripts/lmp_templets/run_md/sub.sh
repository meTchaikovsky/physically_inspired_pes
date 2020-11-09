#!/bin/bash
#SBATCH --job-name=lmp
#SBATCH --partition=TH_HPC2
#SBATCH --nodes=1

mpirun -n 28 lmp < in.input
