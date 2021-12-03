library(tidyverse)
library(vegan)
Phenos = read_tsv("Data/Immunity_Phenotype_list_for_Alex_28_10_2021__TS.txt")
colnames(Phenos)[1] = "ID"
DF = read_rds("Data/Immuno_matrix_postSelection.rds")
DF %>% filter(grepl("32_", ID)) -> DF
sapply(DF$ID, function(x){ str_split(x, "_")[[1]][2] }) -> New_id
DF %>% mutate(ID = New_id) -> DF
Cov = read_tsv("Data/Covariates_LLD&IBD.tsv")
left_join(DF, Cov) %>% mutate(ID = Sample_name) %>% filter(! is.na(Sex)) %>% select(-Sample_name) -> DF
      

Permanova =  function(Distance, Data, Formula){
        adonis2(formula = Formula , data= Data  , na.action = na.exclude, permutations = 2000 ) -> permanova
        as.data.frame(permanova)["Phenotype",] -> permanova
        return( tibble(R2 = permanova$R2, P = permanova$`Pr(>F)`) )

}


DF %>% filter(! is.na(ID)) %>% select_if(~ !any(is.na(.))) -> DF2
All_results = tibble()

print("Testing covariates")

left_join(Phenos, mutate(Cov, ID=Sample_name )) %>% drop_na() %>% arrange(ID) -> Input
DF2 %>% arrange(ID) %>% filter(ID %in% Input$ID) -> DF_i
print(colnames(DF_i))
Distance = vegdist( select(DF_i, -c("ID","Cohort","Probe_n","Sex","Age","plate_id","Disease_status")), method="jaccard" )
#For Sex
Formula = as.formula("Distance ~ plate_id + Age + Sex")
adonis2(formula = Formula , data= Input  , na.action = na.exclude, permutations = 2000 ) -> permanova
as.data.frame(permanova)["Sex",] -> permanova
result = tibble(R2 = permanova$R2[1], P = permanova$`Pr(>F)`[1], Phenotype = "Sex")
#For Age
Formula = as.formula("Distance ~ Sex + plate_id + Age")
adonis2(formula = Formula , data= Input  , na.action = na.exclude, permutations = 2000 ) -> permanova
as.data.frame(permanova)["Age",] -> permanova
result2 = tibble(R2 = permanova$R2[1], P = permanova$`Pr(>F)`[1], Phenotype="Age")
#For plate id
Formula = as.formula("Distance ~ Age+ Sex + plate_id")
adonis2(formula = Formula , data= Input  , na.action = na.exclude, permutations = 2000 ) -> permanova
as.data.frame(permanova)["plate_id",] -> permanova
result3 = tibble(R2 = permanova$R2[1], P = permanova$`Pr(>F)`[1], Phenotype = "plate_id")
rbind(rbind(rbind(result,result2), result3), All_results) -> All_results
print("Starting loop")
for (Pheno in colnames(Phenos)){
        if (Pheno %in% c("ID","antrop_age", "antrop_sex_1_female")){ next }
        Phenos %>% select(c("ID",Pheno)) %>% drop_na() %>% filter(ID %in% DF2$ID) %>% arrange(ID) -> Subset     
        
        DF2  %>% filter(ID %in% Subset$ID) %>% arrange(ID) -> DF25
        DF25  %>% select(-c(ID, Sex, Age, plate_id,Cohort, Probe_n, Disease_status)) -> DF3
        print(colnames(DF3))
        Distance = vegdist( DF3 , method="jaccard")
        
        Subset %>% mutate(Phenotype = as_vector(select(Subset, -ID))) -> Subset
        Subset %>% mutate(Age = DF25$Age, Sex = DF25$Sex, plate_id =DF25$plate_id) -> Subset
        print(Subset)   

        Formula = as.formula("Distance ~ plate_id + Age + Sex + Phenotype")
        Permanova(Distance=Distance, Data=Subset, Formula=Formula) -> result
        result %>% mutate(Phenotype = Pheno) -> result
        rbind(All_results, result) -> All_results
}


write_tsv(All_results, "Results/Permanova.tsv")




