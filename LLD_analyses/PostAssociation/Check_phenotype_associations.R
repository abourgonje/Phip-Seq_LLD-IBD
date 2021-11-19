#S. Andreu-Sanchez
#Check phenotype associations

library(tidyverse)
library(pheatmap)
library(readxl)


###Colors and peptide groups####
#Colors, 25 different colors easily identified by human eye
c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
Groups = tibble(
  Taxa =c(
    "Bacteroides ovatus",                                 
    "Dermatophagoides pteronyssinus (house dust mite)",   
    "Haemophilus influenzae Rd KW20",                     
    "Haemophilus parainfluenzae T3T1",                    
    "Herpes simplex virus",                               
    "Human alphaherpesvirus 2",                           
    "Enterococcus faecalis V583",                         
    "Haemophilus influenzae nontypable strain 3179B",     
    "Lactobacillus fermentum IFO 3956",                   
    "Parabacteroides",                                    
    "Haemophilus influenzae",                             
    "Human herpesvirus 1",                                
    "Lactobacillus acidophilus NCFM",                     
    "Rhinovirus B",                                       
    "Bubalus bubalis (buffalo)",                          
    "Escherichia coli O127:H6 str. E2348/69",             
    "Escherichia coli O157:H7 str. EDL933",               
    "Streptococcus pyogenes",                             
    "Bacteroidales",                                      
    "Bacteroides",                                        
    "Human rhinovirus",                                   
    "Influenza A virus",                                  
    "Streptococcus lutetiensis 033",                      
    "Ovis aries (sheep)",                                 
    "Staphylococcus aureus subsp. aureus MW2",            
    "Streptococcus pneumoniae",                           
    "Enterovirus B",                                      
    "Human alphaherpesvirus 1",                           
    "Staphylococcus aureus",                              
    "Enterovirus",                                        
    "Clostridiales",                                      
    "Staphylococcus aureus subsp. aureus MRSA252",        
    "Bos taurus",                                         
    "Human gammaherpesvirus 4 (EBV)",                     
    "Human betaherpesvirus 5 (CMV)"), 
  Group = c(
    "Bacteroides",                                 
    "Mite)",   
    "Haemophilus",                     
    "Haemophilus",                    
    "Herpesvirus",                               
    "Herpesvirus",                           
    "Enterococcus",                         
    "Haemophilus",     
    "Lactobacillus",                   
    "Parabacteroides",                                    
    "Haemophilus",                             
    "Herpesvirus",                                
    "Lactobacillus",                     
    "Rhinovirus",                                       
    "Mammal",                          
    "Escherichia",             
    "Escherichia",               
    "Streptococcus",                             
    "Bacteroidales",                                      
    "Bacteroides",                                        
    "Rhinovirus",                                   
    "Influenza",                                  
    "Streptococcus",                      
    "Mammal",                                 
    "Staphylococcus",            
    "Streptococcus",                           
    "Enterovirus",                                      
    "Herpesvirus",                           
    "Staphylococcus",                              
    "Enterovirus",                                        
    "Clostridiales",                                      
    "Staphylococcus",        
    "Mammal",                                         
    "Herpesvirus",                     
    "Herpesvirus"                      
  )
)
#######

#Read summary statistics of phenotype associations
ead_tsv("~/Resilio Sync/Antibodies_WIS (1)/Results/Phenotype_associations_Oct31/Pheno_31Oct.txt") -> Antibodies
Antibodies %>% filter(! pheno %in% c("Immune_ctd_strict_greater_1", "Immune_pos_anti_dsDNA", "Immune_pos_CCP", "Immune_pos_SSA_preg", "Immune_ctd_loose_great_0_7") ) %>%
  mutate(FDR = p.adjust(P, "fdr")) %>% arrange(FDR) -> Antibodies
readxl::read_xlsx("~/Resilio Sync/Antibodies_WIS (1)/Data and oligo info/Data_subset/2816-2798-2893_uptodate probes_annotation.xlsx") -> Annotation
Annotation %>% mutate(probe = order ) %>%  select(probe, Taxa, Description) -> Annotation
left_join(Antibodies, Annotation, by="probe") -> Antibodies
#write_tsv(Antibodies, "~/Resilio Sync/Antibodies_WIS (1)/Results/Phenotype_associations_Oct31/Pheno_31Oct_annotated.txt") 

readxl::read_excel("~/Resilio Sync/Antibodies_WIS (1)/Results/Phenotype_associations_Oct31/Pheno_31Oct_annotated_v2_SZ.xlsx",sheet = "Phe_31Oct_annot_v2_fdr0.05") -> Antibodies



#Read annotation made by Sergio, to replace the names of the phenotypes by something more readable and that takes less space
readxl::read_xls("~/Resilio Sync/Antibodies_WIS (1)/Phenotypes/Change_names_heatmap.xls",sheet = 2) -> Replace_n

#Keep only significant associations
Antibodies %>% filter(FDR < 0.05) -> Antibodies_s


#######################################################################
#Get metrix about number or probes and number of phenotypes associated
#Significant associations
Antibodies_s %>% dim()

Antibodies_s %>% group_by(probe) %>% summarise(N = n()) %>% arrange(desc(N)) -> Probe_associations
Antibodies_s %>% group_by(pheno) %>% summarise(N = n()) %>% arrange(desc(N)) -> pheno_associations

#How many associated probes
dim(Probe_associations)[1]
summary(Probe_associations$N)
#How many associated  phenotypes
dim(pheno_associations)[1]

pheno_associations %>% filter(grepl("lymphocytes", pheno))
pheno_associations %>% filter(grepl("neutroph", pheno))
pheno_associations %>% filter(grepl("smk", pheno))
pheno_associations %>% filter(grepl("Allergy", pheno))

#Check smoking
Antibodies_s %>% filter(grepl("smk_now", pheno)) -> SMKNOW
Antibodies_s %>% filter(grepl("smk_full_year", pheno)) -> ever
SMKNOW$probe
length(intersect(SMKNOW$probe,ever$probe))
length(union(SMKNOW$probe,ever$probe))
#Check CCP
Antibodies_s %>% filter(grepl("CCP", pheno))
#Check age and sex
Antibodies_s %>% filter(pheno == "Age") %>% select(-c(se,F,P,FDR,`probe.NA`, `probe.positive`)) %>% group_by(beta>0) %>% summarise(n())
Antibodies_s %>% filter(pheno == "Age") -> Age_associations
write_tsv(Age_associations, "~/Resilio Sync/Antibodies_WIS (1)/Results/Phenotype_associations_Oct31/Age_associations.tsv")
Antibodies_s %>% filter(pheno == "Age") %>% select(-c(se,F,P,FDR,`probe.NA`, `probe.positive`)) %>% print(n=43)
Antibodies_s %>% filter(pheno == "Sex") %>% select(-c(se,F,P,FDR,`probe.NA`, `probe.positive`)) %>% print(n=43)
Antibodies_s %>% filter(pheno == "Sex") -> Sex_associations
write_tsv(Sex_associations, "~/Resilio Sync/Antibodies_WIS (1)/Results/Phenotype_associations_Oct31/Sex_associations.tsv")


###################################################################

#Figure overall barplot

colnames(Replace_n)[1] = "pheno"
left_join(Antibodies_s, Replace_n) -> Phenotypes2


left_join(Phenotypes2, Groups) -> Phenotypes2
Phenotypes2 %>% group_by(category) %>% summarise(N = n()) %>% arrange(N) -> LEVEL

Phenotypes2 %>% mutate(Group = ifelse(is.na(Group), "zOther", Group)) %>%
  mutate(category = fct_rev(factor(category, levels = LEVEL$category ))) -> Phenotype_plot
saveRDS(Phenotype_plot, "~/Input_figures/Fig3A.rds")
Phenotype_plot %>%
ggplot(aes(x=category ,fill=Group)) + geom_bar( ) + theme_bw() + coord_flip()  + scale_fill_manual(values=c25) -> Barplot #scale_y_continuous(trans = "mylog10") 
ggsave("~/Resilio Sync/Antibodies_WIS (1)/Results/Figures_LLD_manuscript/BarPlot_phenotype.pdf",Barplot)


########################################################
Groups = tibble(
  Taxa =c(
    "Streptococcus pneumonia",
    "Streptococcus lutetiensis 033",
    "root",
    "Rhinovirus serotype 2",
    "Rhinovirus B",
    "Pseudomonadaceae",
    "Lactobacillus reuteri DSM 20016",
    "Lactobacillus fermentum IFO 3956",
    "Human Rhinovirus Serotype 2",
    "Human rhinovirus",
    "Human poliovirus",
    "Human gammaherpesvirus 4 (EBV)",
    "Escherichia coli O157:H7 str. EDL933",
    "Escherichia coli O127:H6 str. E2348/69",
    "Enterovirus_Coxsackievirus",
    "Enterovirus B",
    "Enterovirus",
    "Enterococcus faecalis V583",
    "Coxsackievirus A9",
    "Campylobacter phage CP30A",
    "Bos taurus",
    "Bacteroidales",
    "Streptococcus pneumoniae"
  ), Group = c(
    "Streptococcus",
    "Streptococcus",
    "Other",
    "Rhinovirus",
    "Rhinovirus",
    "Pseudomonadaceae",
    "Lactobacillus",
    "Lactobacillus",
    "Rhinovirus",
    "Rhinovirus",
    "Poliovirus",
    "Herpesvirus",
    "Escherichia",
    "Escherichia",
    "Enterovirus",
    "Enterovirus",
    "Enterovirus",
    "Enterococcus",
    "Coxsackievirus",
    "Campylobacter",
    "Bos taurus",
    "Bacteroidales",
    "Streptococcus"),
)

########################################################

#Specific phenotype plots

#Allergy plot
Phenotypes2 %>% filter(category == "Allergies") %>% mutate(Group= ifelse(grepl("Triticum", Taxa), "Plant", Group), Group = ifelse(is.na(Group), "Bacteria", ifelse(Group=="Clostridiales", "Bacteria", Group)) ) -> Fig4C
saveRDS(Fig4C, "~/Input_figures/Fig4C.rds")

Fig4C %>%
  ggplot(aes(x=pheno ,fill=Group)) + geom_bar( ) + theme_bw() + coord_flip()  + scale_fill_manual(values=c25)  

#auto-abs plot
Phenotypes2 %>% filter(category == "Autoimmune") %>% 
  mutate(Group = ifelse(! is.na(Group), Group, ifelse(grepl("allergen", probe.protein), "Mammal", ifelse(grepl("agilent", probe), "Bacteria", "Virus" ) )  )  ) -> Phenos_auto
  Phenos_auto %>% mutate( Group = ifelse(Group %in% c("Bacteroidales", "Bacteroides", "Haemophilus","Lactobacillus","Parabacteroides", "Staphylococcus", "Streptococcus"), "Bacteria", ifelse(Group == "Herpesvirus","Virus",Group)  ) ) -> Phenos_auto
  saveRDS(Phenos_auto, "~/Input_figures/Fig4D.rds")
  
  Phenos_auto %>% ggplot(aes(x=pheno ,fill=Group)) + geom_bar( ) + theme_bw() + coord_flip()  + scale_fill_manual(values=c25)  



##Smoking-specific plots
#Transfer is a matrix of smoking probes
read_tsv("~/Desktop/Transfer") -> AB2
info = read_csv("~/Resilio Sync/Antibodies_WIS (1)/Data and oligo info/Data_subset/2821probes_annotation.csv")
info %>% mutate(probe = order) %>% select(probe, Taxa, Description) -> info
AB2 %>% drop_na() -> AB2
#Divide in two groups
AB2 %>% filter(smk_now == 2) %>% select(-smk_now) -> Smoker
AB2 %>% filter(smk_now == 1) %>% select(-smk_now) -> Non_Smoker
#Get prevalences
apply(Smoker,2,function(x){ x = x[! is.na(x)] ;  sum(x)/length(x) } )  %>%  as_vector() -> Prevalence_smoker
apply(Non_Smoker,2, function(x){ x = x[! is.na(x)] ;  sum(x)/length(x) } ) %>% as_vector() -> Prevalence_nonsmoker
#Combine in a wide format
tibble(Prevalence = Prevalence_smoker, Smoker="Yes") -> smk_t
tibble(Prevalence = Prevalence_nonsmoker, Smoker="No") -> nsmk_t
cbind(smk_t, nsmk_t) -> To_plot
#Plot 3, smoking prevalence comparison
colnames(To_plot)= c("Probe prevalence\nSmokers", "R1", "Probe prevalence\nnon-smokers", "R2")
select(To_plot, -c("R1", "R2")) %>% mutate(probe = colnames(select(AB2, -smk_now)) ) %>% as_tibble() -> To_plot
left_join(left_join(To_plot, Annotation2), Groups) -> To_plot
To_plot %>% ggplot(aes(x=`Probe prevalence\nSmokers`, y=`Probe prevalence\nnon-smokers`, col=Group)) + theme_bw() + geom_point() + geom_abline() + scale_color_manual(values=c25) -> Smoke_plot1
ggsave(Smoke_plot1 ,"~/Resilio Sync/Antibodies_WIS (1)/Results/Figures_LLD_manuscript/Smoking_enriched_probes.pdf")
#Plot 4. Heatmap of smoking probes
AB2 %>% select(-smk_now) %>% as.data.frame()-> To_plot
AB2 %>% mutate(Smoker = ifelse(smk_now == 1, "No", "Yes") ) %>% select(Smoker) %>% as.data.frame() -> Annotation
rownames(To_plot) = rownames(Annotation)
pheatmap::pheatmap(To_plot, annotation_row = Annotation, show_rownames = F)
#Add annotation
readxl::read_xlsx("~/Resilio Sync/Antibodies_WIS (1)/Results/PhenoAllergMerged.Oct4_SZ.xlsx",sheet = "FDR0.05") %>%
  filter(probe %in% colnames(AB2), pheno == "smk_now") %>% select(probe, Taxa) -> Annotation2
left_join(Annotation2, Groups) -> Annotation2
Annotation2 %>% select(-Taxa) %>% as.data.frame() %>% column_to_rownames("probe") -> Annotation2
Colors = list(Smoker = c("Yes"= "Red", "No" = "Blue") )
pheatmap::pheatmap(To_plot, annotation_row = Annotation,annotation_col = Annotation2, show_rownames = F, annotation_colors = Colors )

#Plot distances between smoking probes
Distance_smoke = "~/Desktop/Smoking_probes.mat"
readxl::read_xlsx("~/Resilio Sync/Antibodies_WIS (1)/Results/PhenoAllergMerged.Oct4_SZ.xlsx",sheet = "FDR0.05") %>%
  filter(probe %in% colnames(AB2), pheno == "smk_now") %>% select(probe, Taxa) -> Annotation2
left_join(Annotation2, Groups) -> Annotation2

read.table(Distance_smoke, row.names=1, skip=1) -> Distance
cbind(smk_t, nsmk_t) -> To_plot2 ; colnames(To_plot2)= c("Probe prevalence\nSmokers", "R1", "Probe prevalence\nnon-smokers", "R2") ;  select(To_plot2, -c("R1", "R2")) %>% mutate(probe = colnames(select(AB2, -smk_now)) ) %>% as_tibble() -> To_plot2
left_join(left_join(To_plot2, Annotation2), Groups) -> To_plot2

lapply(rownames(Distance), function(x){ str_split(x,":")[[1]][1] } ) %>% unlist() -> Names
rownames(Distance) = paste(To_plot2$Group, seq(1, length(Names)))
colnames(Distance) = rownames(Distance)
pheatmap::pheatmap(Distance, color=viridis::viridis(10))



