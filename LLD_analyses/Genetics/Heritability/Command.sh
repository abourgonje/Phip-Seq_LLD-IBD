gcta="/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/Programs/gcta_1.93.2beta/gcta64"
Matrix="/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/1.Relationship_matrix/GRM_LLD"

Pheno=$1
Pheno_n=$2


mkdir /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/2.GREML/H/$Pheno_n
Out=/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/2.GREML/H/$Pheno_n/Output



$gcta  --grm $Matrix --pheno $Pheno --reml --out $Out --thread-num 1 --covar  /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/2.GREML/Categorical_covariates.covar --qcovar /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/2.GREML/Numerical_covariates.covar

