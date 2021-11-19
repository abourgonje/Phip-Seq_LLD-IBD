#S. Andreu-Sanchez
#Explore HLA associations

library(tidyverse)
library(readxl)
library(patchwork)
c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
         "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

#Read data and annotation
read_tsv("~/Resilio Sync/Antibodies_WIS (1)/Results/HLA_meta/MetaHLA_Significant_FDR.tsv") -> df
Annotation = readxl::read_excel("~/Resilio Sync/Antibodies_WIS (1)/Data and oligo info/Data_subset/2816-2798-2893_uptodate probes_annotation.xlsx") %>% filter(order %in% df$probe) 
Annotation %>% select(order, aa_seq, Taxa, Description, Assoc_pheno_top) -> Annotation
colnames(Annotation)[1] = "probe"
#General checks
df %>% group_by(Type) %>% summarise(n())
df %>% filter(Bonferroni_ind < 0.05) %>% group_by(Type) %>% summarise(n())
df %>% filter(Bonferroni_ind < 0.05) -> df
df %>% group_by(probe) %>% summarise(n())
df %>% select( - c("Unnamed: 0","IBD_A1", "IBD_A2","IBD_SE", "LLD_A1", "LLD_A2", "LLD_SE"   )) -> df
###Odd min value
summary(df$IBD_BETA)
summary(df$LLD_BETA)
summary(df$OR) #This looks normal


df %>% group_by(Type) %>% summarise(N= n()) %>% ggplot(aes(x="", y=N, fill=Type)) + geom_bar(stat = "identity",width=1, color="white") +
  coord_polar("y", start=0) + theme_void() 

#Work only on HLA alleles instead of SNP/aa, add some new columns with info about HLA I/II, chain and allele
df %>% filter(Type == "HLA") -> df_genes
df_genes %>% mutate(HLA_class = ifelse(grepl("D", SNP), "II", "I"),
                    Gene_class = ifelse(grepl("_A_", SNP), "A", ifelse(grepl("_B_", SNP), "B",ifelse(grepl("_C_", SNP), "C", "other"))) ) %>%
mutate(Gene_class = ifelse( grepl("DQ", SNP), "DQ", ifelse(grepl("DR", SNP), "DR", ifelse(grepl("DP", SNP), "DP",  Gene_class)))) %>%
  mutate(Chain = ifelse(HLA_class == "I", "Alpha", ifelse(grepl("B1", SNP), "Beta", "Alpha"))  ) -> df_genes
df_genes %>% select(-c(N,P, `P(R)`,`OR(R)`, Q, I, IBD_BETA, LLD_BETA )) -> df_genes
df_genes %>% arrange(FDR) -> df_genes
df_genes %>% group_by(HLA_class) %>% summarise(n())

#General plots
df_genes %>% ggplot(aes(x=log10(OR), y=-log10(FDR), shape=Chain, col=Gene_class, fill=Chain)) + geom_point() + theme_bw() +
  facet_wrap(~HLA_class) + scale_color_manual(values=c25)
df_genes %>% ggplot(aes(x="", y=, fill=HLA_class)) + geom_bar(width=1, color="white") + coord_polar("y", start=0) + theme_void() -> HLA_class_plot

#Seems that most associations are in beta chain
df_genes %>% filter(HLA_class == "II") %>% group_by(Chain) %>% summarise(n())
df_genes %>% filter(HLA_class == "II") %>% ggplot(aes(x="", y=, fill=Chain)) + geom_bar(width=1, color="white") + coord_polar("y", start=0) + theme_void() -> Chain_plot
#Most associations are in DQ and DR, no DP
df_genes %>% filter(HLA_class == "II") %>% group_by(Chain, Gene_class) %>% summarise(N = n()) -> Counts_chain_gene
Counts_chain_gene %>% mutate(N = ifelse(Chain=="Alpha", N/270, N/800 ) ) %>%
ggplot(aes(x="", y=N, fill=Gene_class)) + geom_bar(stat = "identity",width=1, color="white") +
 facet_wrap(~Chain)  +coord_polar("y", start=0) + theme_void() -> Gene_class_plot
#DR has plenty assocaitions in beta chains, but none in alpha chains
HLA_class_plot + Chain_plot / Gene_class_plot



df_genes %>% group_by(probe, HLA_class) %>% summarise(N = n() ) %>% arrange(desc(N)) #338
df_genes %>% group_by(probe) %>% summarise(N = n() )  #300
#Meaning that there are probes that associated both with HLA I and II
df_genes %>% group_by(probe, HLA_class, OR>1) %>% summarise(N = n() ) %>% arrange(desc(N)) #1,157 


df_genes %>% filter(probe=="agilent_235068")  %>% ggplot(aes(x = SNP, y =OR, col=Chain )) +  geom_point(stat="identity") + 
  facet_wrap(~Gene_class, scales = "free") + theme_bw() + coord_flip() + ggtitle("agilent_235068")
df_genes %>% filter(probe=="twist_88894")  %>% ggplot(aes(x = SNP, y =OR, col=Chain )) + coord_flip() + 
    geom_point()  + #geom_text(aes(label = round(OR,2), hjust = -0.5)) + 
  facet_wrap(~Gene_class, scales = "free") + theme_bw() + theme(panel.grid.major.y = element_blank())  + ggtitle("twist_88894")
df_genes %>% filter(probe=="agilent_241118")  %>% ggplot(aes(x = SNP, y =OR, col=Chain )) + coord_flip() +  geom_point()  +
  facet_wrap(~Gene_class, scales = "free") + theme_bw() + theme(panel.grid.major.y = element_blank())  + ggtitle("agilent_241118")
df_genes %>% filter(probe=="agilent_9680")  %>% ggplot(aes(x = SNP, y =OR, col=Chain )) + coord_flip() +  geom_point()  +
  facet_wrap(~Gene_class, scales = "free") + theme_bw() + theme(panel.grid.major.y = element_blank())  + ggtitle("agilent_9680")

df %>% filter(Type == "HLA" ) %>% filter(probe == "twist_88894") %>% arrange(P) %>% select(SNP, P, OR)

##Check annotation
Annotation = readxl::read_excel("~/Resilio Sync/Antibodies_WIS (1)/Data and oligo info/Data_subset/2798-2893_probes_annotation.xlsx")
Annotation %>% mutate(probe = order) %>% select(probe, Taxa, Description,aa_seq) -> Annotation
left_join(df_genes, Annotation) -> df_genes_annotated

df_genes_annotated %>% filter(HLA_class == "II", Gene_class == "DR" ) %>% group_by(SNP) %>% summarise(n())
df_genes_annotated %>% filter(HLA_class == "I", Gene_class == "A" ) %>% group_by(SNP) %>% summarise(n())
df_genes_annotated %>% filter(HLA_class == "I", Gene_class == "B" ) %>% group_by(SNP) %>% summarise(n())
df_genes_annotated %>% filter(HLA_class == "I", Gene_class == "C" ) %>% group_by(SNP) %>% summarise(n())
df_genes_annotated %>% filter(HLA_class == "II", Gene_class == "DQ", Chain =="Alpha" ) %>% group_by(SNP) %>% summarise(n())
df_genes_annotated %>% filter(HLA_class == "II", Gene_class == "DQ", Chain =="Beta" ) %>% group_by(SNP) %>% summarise(n())
##


#####Compare predicted HLA-peptide affinities with associations


df_genes_annotated %>% select(probe, aa_seq) %>% distinct() -> Peptides 
write_tsv(Peptides, "~/Desktop/Peptides.fa")



#Check affinities
read_csv("~/Downloads/All_results.csv") -> Affinities_DRB1
df_genes %>% filter(Gene_class == "DR", HLA_class == "II") -> To_merge
sapply(To_merge$SNP, function(x){ str_replace(x, "HLA_", "") } ) -> Variants

To_merge$SNP = Variants
##Clumping of top variants
To_merge %>% filter(HLA_class == "II", Chain == "Beta") %>% select(-c("CHR","BP","A1","A2","OR", "SNP")) %>% group_by(probe, Gene_class) %>% summarise(n()) %>% arrange(probe) %>% summarise(n()) -> Associated_to_both
To_merge %>% filter(probe %in% Associated_to_both$probe)
New_ones = tibble()
for (P in Associated_to_both$probe){
  To_merge %>% filter(HLA_class == "II", Chain == "Beta") %>% filter(probe == P) %>% arrange(FDR) %>% head(n=1) -> C
  rbind(New_ones, C) -> New_ones
}
To_merge  %>% filter(! probe %in% Associated_to_both$probe) -> Rest
rbind(Rest, New_ones) -> Clumped
Clumped %>% filter(HLA_class == "II", Chain == "Beta") -> Clumped
Clumped %>% mutate(ID = paste(SNP, probe) ) -> Clumped
##

To_merge %>% mutate(MHC = Variants, Identity=probe ) %>% select(-c(SNP, CHR, BP, A1, A2, probe) )  -> To_merge
left_join(Affinities_DRB1, To_merge, by=c("MHC", "Identity") ) -> Merged_v
Merged_v %>% mutate(ID = paste(MHC ,Identity)) %>%  filter(! ID %in% Clumped$ID) -> Merged_v
Merged_v %>% filter(is.na(OR)) -> Not_associated
Merged_v %>% filter(!is.na(OR)) -> Associated


Merged_v %>% mutate(HLA_asociation = ifelse(X1 %in% Not_associated$X1, "Not_significant", "Significant") ) -> Merged_v

Merged_v %>% mutate(Binding = ifelse(BindLevel == "<=NB" , "no", "yes")) %>% group_by(Binding, HLA_asociation) %>% summarise(N = n()) -> R_n
R_n %>% spread(Binding, N) %>% as.data.frame() %>% column_to_rownames("HLA_asociation") -> test_table
fisher.test(test_table)
Merged_v %>% mutate(Binding = ifelse(BindLevel == "<=SB" , "yes", "no")) %>% group_by(Binding, HLA_asociation) %>% summarise(N = n()) -> R_s
R_s %>% spread(Binding, N) %>% as.data.frame() %>% column_to_rownames("HLA_asociation") -> test_table_s
fisher.test(test_table_s)

Merged_v %>% mutate(Binding = ifelse(BindLevel == "<=SB" , "yes", "no")) %>% ggplot(aes(x=HLA_asociation, fill=BindLevel)) + 
  geom_bar(position="fill")
Merged_v %>% mutate(Strong_binding = ifelse(BindLevel == "<=SB" , "yes", "no")) %>% ggplot(aes(x=HLA_asociation, fill=Strong_binding)) + 
  geom_bar(position="fill")
