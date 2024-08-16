#!/bin/bash
#SBATCH --job-name=ljfit
#SBATCH --partition=infinite
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1300

python Script_fit.py

# We succeeded, reset trap and clean up normally.
trap - EXIT
exit 0
