#S. Andreu-Sanchez
#Given a matrix of 0/1s and genetic correlations, it compares the correlations of the 0/1 profiles with the genetic correlations

library(tidyverse)
library(readxl)
library(corrplot)
library(vegan)

set.seed(99)


#Annotation
readxl::read_excel("Resilio Sync/Antibodies_WIS (1)/Data and oligo info/Data_subset/2798-2893_probes_annotation.xlsx") -> Annotation
#Heritability
read_tsv("Resilio Sync/Antibodies_WIS (1)/Results/Genetic_correlation/GeneticCorrelation_LLD_above50heritability.tsv", col_names = F) -> GC
colnames(GC) = c("Probe", "Probe2", "Genetic_correlation", "SE" )


#Make a correlation matrix of r_g
GC %>% mutate(Genetic_correlation = ifelse(Probe == Probe2, 1, Genetic_correlation)) %>% mutate(Genetic_correlation = ifelse(SE > Genetic_correlation, 0, Genetic_correlation) ) -> GC
GC %>% select(Probe, Probe2, Genetic_correlation) %>% spread(Probe2, Genetic_correlation) %>% as.data.frame() %>% column_to_rownames("Probe") -> GC_matrix
GC_matrix[is.na(GC_matrix)] = 1
#Make 6 clusters using KMeans
kmeans(GC_matrix, 6)$cluster %>% as.factor() -> clusters
data.frame(Probe = names(clusters), Cluster=as.vector(clusters)) -> clusters2
#Save clusters
pdf("~/Resilio Sync/Antibodies_WIS (1)/Results/Genetic_correlation/Heatmap_genetic_Correlation.pdf",   width = 4, height = 4)
pheatmap::pheatmap(GC_matrix, annotation_row = column_to_rownames(clusters2, "Probe"))
dev.off()

#Pull description of each cluster
left_join(GC, clusters2) %>% filter(! Cluster %in% c(5,4,1)) -> Keep_GC
colnames(clusters2) = c("Probe2", "Cluster2")
left_join(Keep_GC, clusters2 ) %>% filter(Cluster == Cluster2) -> Keep_GC
for (C in unique(Keep_GC$Cluster)){
  print(paste(c("Cluster : ", as.character(C)), collapse="" ))
  Keep_GC %>% filter(Cluster == C) -> CHECK
  Annotation %>% filter(order %in% c(CHECK$Probe, CHECK$Probe2)) -> CHECK2
  CHECK2 %>% select(order,Taxa, Description) %>% print()
}
#3: Plasmodium (3 probes)
#2: Clostridiales, Lactobacillus phage, S. aureus and wheat 
#2: Streptococcus, aureus, pneuomoniae, pyogenes
### Plasmodium is weird, dig
clusters2 %>% filter(Cluster2 == 1)
Annotation %>% filter(order %in% filter(clusters2,Cluster2 == 1)$Probe2 ) %>% select(order, Prevalence_LLD, Taxa, Description, aa_seq)
Annotation %>% filter(order %in% filter(clusters2,Cluster2 == 6)$Probe2 ) %>% select(order, Prevalence_LLD, Taxa, Description, aa_seq)

########################################################
#Compare co-abundance with Genetic_correlation##########
########################################################

Data = readRDS("~/Desktop/Immuno_matrix_postSelection.rds") ;  process( filter(Data, grepl("32_", ID))) -> Data_LLD
Data_LLD %>% select( unique(GC$Probe, GC$Probe2) ) -> Subset_LLD
cor(Subset_LLD) -> Correlations_LLD

rownames(Correlations_LLD)
rownames(GC_matrix)
test = hclust(dist(GC_matrix))
ORDER = rownames(GC_matrix)[test$order]
Correlations_LLD[ORDER,ORDER] ->  Correlations_LLD
GC_matrix[ORDER,ORDER] ->  GC_matrix
New_matrix = Correlations_LLD
New_matrix[lower.tri(Correlations_LLD)] <- as.matrix(GC_matrix)[lower.tri(as.matrix(GC_matrix))]

corrplot(New_matrix, diag=FALSE, tl.col="black") -> Plot_m
corrMatOrder(corr = Plot_m[test$order], order = names(GC_matrix)[test$order]  )

#mantel test
mantel( as.dist(Correlations_LLD) , as.dist(GC_matrix), method = "spearman", permutations = 9999, na.rm = TRUE)








