#!/bin/sh
#SBATCH  --job-name=Family.job
#SBATCH --output=logs/family.out
#SBATCH --error=logs/family.err
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1

ml RPlus/4.0.3-foss-2018b-v21.08.1
Rscript Consistency/Family_analysis.R

