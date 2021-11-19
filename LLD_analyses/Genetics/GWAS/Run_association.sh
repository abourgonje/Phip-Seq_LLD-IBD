#!/bin/sh
#SBATCH  --job-name=GWAS.job
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --nodes=1


ml PLINK/1.9-beta6-20190617



Name=$1
Out=$2
Name_p=$3

if [ "$Name" = "Number_probes" ]; then
plink --bfile ../Genetic_Preprocess/2_Relatedness_ethnicity/LLD  --covar Covariates.txt --linear --out $Out/$Name --fam $Name_p

else
plink --bfile ../Genetic_Preprocess/2_Relatedness_ethnicity/LLD  --covar Covariates.txt --logistic --out $Out/$Name --fam $Name_p #--reference-allele #--hide-covar

fi
