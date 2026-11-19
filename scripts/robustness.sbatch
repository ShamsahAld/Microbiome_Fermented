#!/bin/bash
#SBATCH --job-name=robustness
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/robustness_%j.log
cd "$SLURM_SUBMIT_DIR"
pixi run python scripts/robustness.py
