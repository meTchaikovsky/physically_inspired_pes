#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=debug2
#SBATCH --nodes=1

python train_model.py
