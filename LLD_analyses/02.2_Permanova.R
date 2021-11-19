library(tidyverse)
library(vegan)

Phenos = read_tsv("") #Path to Phenotypes
colnames(Phenos)[1] = "ID"
DF = read_rds("Data/Immuno_matrix_postSelection.rds") #Path to the matrix generated in 01
DF %>% filter(grepl("32_", ID)) -> DF #Selection of LLD cohort
sapply(DF$ID, function(x){ str_split(x, "_")[[1]][2] }) -> New_id
DF %>% mutate(ID = New_id) -> DF #Remove the prefix
Cov = read_tsv("") #Path to covariates
left_join(DF, select(Cov, c(ID, Sample_name))) %>% mutate(ID = Sample_name) %>% select(-Sample_name) -> DF
      

Permanova =  function(Distance, Data, Formula){
	adonis2(formula = Formula , data= Data  , na.action = na.exclude, permutations = 2000 ) -> permanova
	return( tibble(R2 = permanova$R2[1], P = permanova$`Pr(>F)`[1]) )

}


DF %>% filter(! is.na(ID)) %>% select_if(~ !any(is.na(.))) -> DF2
All_results = tibble()
for (Pheno in colnames(Phenos)){
	if (Pheno == "ID"){ next }
	print(Pheno)
	Phenos %>% select(c("ID",Pheno)) %>% drop_na() %>% filter(ID %in% DF2$ID) %>% arrange(ID) -> Subset	

	DF2 %>% filter(ID %in% Subset$ID) %>% arrange(ID) %>% select(-ID) -> DF3
	Distance = vegdist( DF3 , method="jaccard")

	Subset %>% mutate(Phenotype = as_vector(select(Subset, -ID))) -> Subset

	Formula = as.formula("Distance ~ Phenotype")
	Permanova(Distance=Distance, Data=Subset, Formula=Formula) -> result
	result %>% mutate(Phenotype = Pheno) -> result
	rbind(All_results, result) -> All_results
}

write_tsv(All_results, "Results/Permanova.tsv")





