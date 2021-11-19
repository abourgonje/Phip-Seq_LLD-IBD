import sys
from subprocess import call
from pathlib import Path

counter = 0 
batch = []
for path in Path("Phenotypes").glob("*"):
	Name = path.name
	batch.append(str(path))
	counter  +=  1
	if counter == 500:
		script = "Batch_{C}_temp".format(C=str(counter))
		with open(script, "w") as F:
			F.write("""#!/bin/sh
#SBATCH  --job-name=Gencorr.job
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=5G
#SBATCH --nodes=1
gcta="/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/Programs/gcta_1.93.2beta/gcta64"
Matrix="/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/1.Relationship_matrix/GRM_LLD"\n
""")
			for pheno_name in batch:
				Name = pheno_name.split("/")[-1]
				F.write("$gcta  --reml-bivar --grm $Matrix  --pheno {Pheno}  --out Results/{N} --covar  /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/2.GREML/Categorical_covariates.covar --qcovar /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/2.GREML/Numerical_covariates.covar\n".format(Pheno=pheno_name, N=Name))
		exit(script)
		call("sbatch "+ script, shell=True)	
		batch = []
		call("rm "+ script, shell=True)			




