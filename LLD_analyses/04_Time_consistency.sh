#!/bin/sh
#SBATCH  --job-name=Time_consitency.job
#SBATCH --output=logs/Time_consistency.out
#SBATCH --error=logs/Time_consitency.err
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --nodes=1


ml RPlus/4.0.3-foss-2018b-v21.08.1
Rscript Consistency/Consitency_analysis.R
