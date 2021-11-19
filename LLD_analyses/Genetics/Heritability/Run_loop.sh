#!/bin/sh
#SBATCH  --job-name=heritability.job
#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --nodes=1


for F in /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/2.GREML/Phenotypes/* ; do
	n=$(basename $F)
	bash Command.sh $F $n
done 


