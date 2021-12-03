#S. Andreu-Sanchez
#Uses WGCNA in a presence/absence profile to identify modules of highly correlated peptides
#Also: Module correlation. Module composition plotting. Network plotting.
#Sequence similarity clustring

library(WGCNA) # https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/OverviewWGCNA.pdf
library(tidyverse)
library(readxl)
library(igraph)

set.seed(1299)

#Colors, 25 different colors easily identified by human eye
c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "gray70", "khaki2",
        "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
#Functions for WGCNA part
process =  function(Data){
        #cleanup of the data before other functions can be applied
        Data %>% select(-ID) -> Data
        Data[ , apply(Data, 2, function(x) !any(is.na(x)))] -> Data
        return(Data)
}
choose_power = function(Data){
        #Choosing ideal power function.
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold(Data, powerVector = powers, verbose = 5, corFnc = cor )

        sizeGrWindow(9, 5)
        par(mfrow = c(1,2));
        cex1 = 0.9;

        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
                main = paste("Scale independence"));
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        labels=powers,cex=cex1,col="red");
        abline(h=0.90,col="red")

        plot(sft$fitIndices[,1], sft$fitIndices[,5],
                xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
                main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}
Do_analysis = function(Data, k=7, cutoff=30){
        #Build network in one command. 
        net = blockwiseModules(Data, power = k,corType = "pearson",
                       TOMType = "signed", minModuleSize = cutoff,
                       reassignThreshold = 0, mergeCutHeight = 0.5,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = T, randomSeed = 1234,
                       verbose = 3)

        #Distance matrix is saved in BlockwiseTOM-block.1.RData ; can be loaded by load(lockwiseTOM-block.1.RData)
        load("BlockwiseTOM-block.1.RData") #It is called TOM

        #Data of interest
        moduleLabels = net$colors  #classification in labels
        moduleColors = labels2colors(net$colors) #color for the labels in figure
        MEs = net$MEs #Eigengenes
        geneTree = net$dendrograms[[1]] #dendogram
        #Overview of clusters
        table(moduleLabels) -> overview
        #Plot tree
        plotDendroAndColors(geneTree, moduleColors, "Module colors" ,dendroLabels = F, addGuide = T) -> Plot_tree
        #Visualize similarities between the eigen representation of the modules
        Correlation_eigen =  cor(MEs)
        Correlation_eigen[Correlation_eigen>0.999] = NA
        pheatmap::pheatmap( Correlation_eigen ) -> Plot_heatmap
        #Plotting distance between all probes,and coloring by module
        igraph::graph.adjacency(TOM) -> graph_object
        plotTOM = as.matrix(TOM)
        diag(plotTOM) = NA;
        sizeGrWindow(9,9)
        TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all probes") -> Plot_comparison
        tibble(Probe = colnames(Data), Cluster = moduleLabels) -> Probes_and_cluster
        
        return(list(moduleLabels, moduleColors, MEs, geneTree, Plot_tree, Plot_comparison,Probes_and_cluster, net))
}       


#Read data and remove probes with NA
Data = readRDS("~/Desktop/Immuno_matrix_postSelection.rds") #Complete matrix of 0/1s + IDs, script 01
process(Data) -> Data_all
process( filter(Data, grepl("32_", ID))) -> Data_LLD #Removal of samples that are not LLD. Might consider removing longitudinal samples.
read_csv("~/Resilio Sync/Antibodies_WIS (1)/Data and oligo info/Data_subset//2821probes_annotation.csv") -> Annotation #Peptide annotation
colnames(Annotation)[1] = "Probe"

###OVERALL (LLD+IBD)###
#Check power
choose_power(Data_all) #It is recommended to use a number aboce R2 of 0.9, the highest possible before saturation, which in this case is 7
Do_analysis(Data_all, 7) -> results_all
Probes_and_cluster = results_all[[7]]

#Check which are the probes in the same groups
Probes_and_cluster %>% filter(Cluster != 0) -> Groups
left_join(Groups, Annotation) -> Groups
#1 --> CMV , 2 ---> Flagellin, 3 --> Bacteria (allergen?)
Groups %>% select(Cluster, Taxa, Description,`comments, synonim, blast, details`) %>% filter(Cluster == 1) %>% print(n=100) #Cluster CMV, turquoise
Groups %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 2) %>% print(n=100) #Cluster CMV, turquoise
Groups %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 3) %>% print(n=100) #Cluster CMV, turquoise

##LLD##
choose_power(Data_LLD)
Do_analysis(Data_LLD, 7) -> results_LLD
Probes_and_cluster_LLD = results_LLD[[7]]
Probes_and_cluster_LLD %>% filter(Cluster != 0) -> Groups_LLD ;left_join(Groups_LLD, Annotation) -> Groups_LLD

Groups_LLD %>% select(Cluster, Taxa, Description,`comments, synonim, blast, details`) %>% filter(Cluster == 1) %>% print(n=100) #Cluster CMV, turquoise
Groups_LLD %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 2) %>% print(n=100) #Cluster CMV, turquoise
Groups_LLD %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 3) %>% print(n=100) #Cluster CMV, turquoise

results_LLD[[2]][results_LLD[[1]] == 1] #light blue
results_LLD[[2]][results_LLD[[1]] == 2] #dark blue
results_LLD[[2]][results_LLD[[1]] == 3] #brown\

Groups_LLD %>% select(Cluster, Probe,aa_seq) %>% filter(Cluster ==3) -> Cluster_heterogenous
Cluster_heterogenous %>% write_tsv("Desktop/Cluster_heterogenous.tsv")
write_tsv(Groups_LLD, path = "~/Resilio Sync/Antibodies_WIS (1)/Results/Network_analysis/Module_belonging_LLD.tsv")

#Repeat on different cluster number
#With 20 there is one additional module
Do_analysis(Data_LLD, 7, 20)  -> results_LLD2
Probes_and_cluster_LLD2 = results_LLD2[[7]] ; Probes_and_cluster_LLD2 %>% filter(Cluster != 0) -> Groups_LLD2 ;left_join(Groups_LLD2, Annotation) -> Groups_LLD2
Groups_LLD2 %>% select(Cluster, Taxa, Description, `comments, synonim, blast, details`) %>% filter(Cluster == 4) %>% print(n=100)
#New cluster consists of 26 probes. No clear taxa or domain, mainly bugs
results_LLD2[[2]][results_LLD2[[1]] == 4] #yellow cluster

##################Analysis Paper ####################33

##1. Identify modules with at least 10 probes per module
#With 10
Do_analysis(Data_LLD, 7, 10)  -> results_LLD3
Net =  results_LLD3[[8]]
Probes_and_cluster_LLD3 = results_LLD3[[7]] ; Probes_and_cluster_LLD3 %>% filter(Cluster != 0) -> Groups_LLD3 ;left_join(Groups_LLD3, Annotation) -> Groups_LLD3
table(results_LLD3[[2]]) #22 clusters
#Save
write_tsv(Groups_LLD3, "Cluster_AtLeast10.tsv")


#############Correlation between modules using eigengenes
Net$MEs %>% as_tibble() -> Eigengenes
cor(Eigengenes) -> Cor_modules
Cor_eig = tibble()
for (Ei in colnames(Eigengenes)){
        for (Ei2 in colnames(Eigengenes)){
                if (Ei == Ei2){ next }
                ID =paste(sort(c(Ei, Ei2)), collapse="-")
                if ( ID %in% Cor_eig$ID ){ next }
                cor.test(as_vector(select(Eigengenes, Ei)), as_vector(select(Eigengenes, Ei2))) -> test
                rbind(Cor_eig, tibble(ID=ID, P=test$p.value, Estimate=test$estimate)) -> Cor_eig 
                
        }
}
Cor_eig %>% mutate(P_b =p.adjust(P, "bonferroni")) %>%arrange(P_b) %>% filter(P_b < 0.05 ) %>% filter(!grepl("ME0", ID))
cor.test(Eigengenes$ME18, Eigengenes$ME17) ; cor.test(Eigengenes$ME18, Eigengenes$ME14)
Cor_eig %>% mutate(P_b =p.adjust(P, "bonferroni")) %>%arrange(P_b) %>% filter(! grepl("ME0", ID)) %>% filter(P_b < 0.05)

pheatmap::pheatmap(Cor_modules) #Modules 18, 17,14 and ?0? are related
#######################################

#####################################################
###Plot piechart of composition per cluster##########
#####################################################
#Load cluster with annotation
readxl::read_xlsx("Cluster_AtLeast10.xlsx") -> Annotation_groups
Annotation_groups %>% group_by(Cluster, High_taxonomy) %>% summarise(N = n()) -> Composition_cluster
New_cluster = tibble()
Colors_do = c25[2: length(c25)] ; Colors_do[6] = c25[11]
#Make fractions of each category per cluster
for (C in unique(Composition_cluster$Cluster)){
        Composition_cluster %>% filter(Cluster == C) -> subset_composition
        sum(subset_composition$N) -> Total 
        subset_composition %>% mutate(Fraction = N/Total ) -> subset_composition
        rbind(New_cluster, subset_composition) -> New_cluster
        
}
#Plot
ggplot(New_cluster, aes(x="", y=Fraction, fill=High_taxonomy)) + geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) + theme_void() + facet_wrap(~Cluster) + scale_fill_manual(values=Colors_do) -> Piecharts
ggsave("Piecharts_moduleComposition.pdf")

##################################
##Network visualization###########
##################################

Data_LLD %>% select(Groups_LLD3$Probe) -> Module_info
Network_labels = tibble(Name = names(Net$colors), Cluster = Net$colors,  Color =labels2colors(Net$colors) )
Network_labels %>% filter(Name %in% colnames(Module_info)) -> Network_labels

dev.off()


g <- graph.adjacency(
        as.matrix(as.dist(cor(Module_info, method="pearson"))),
        mode="undirected",
        weighted=TRUE,
        diag=FALSE
)
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
E(g)[which(E(g)$weight<0)]$color <- "darkblue"
E(g)[which(E(g)$weight>0)]$color <- "darkred"
E(g)$weight <- abs(E(g)$weight)
g <- delete_edges(g, E(g)[which(E(g)$weight<0.3)])
V(g)$name <- Network_labels$Cluster
V(g)$shape <- "circle" ; V(g)$color <- "skyblue" ; V(g)$vertex.frame.color <- "white"
V(g)$label.cex = Network_labels$Cluster

scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(Module_info, 1, mean)) ) * 10 #Scaled by prevalnece
edgeweights <- E(g)$weight * 2.0
mst <- mst(g, algorithm="prim")

mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=as.vector(Network_labels$Cluster ))

#Color choosing so that they are consistent with the ones used for piecharts
left_join(tibble(Probe = colnames(Module_info)) ,  select(Annotation_groups, c(Probe, High_taxonomy))) -> ForColors
as.factor(ForColors$High_taxonomy) -> ForColors2
Colors_do[1: length(levels(ForColors2))] -> Colors_do
tibble( High_taxonomy = levels(ForColors2), Color = Colors_do ) -> Colors
left_join(ForColors , Colors) -> Colors
V(mst)$color <- Colors$Color
#####
plot(mst,
     layout=layout.fruchterman.reingold,
     edge.curved=TRUE,
     vertex.size=vSizes,
     vertex.label.dist=0,
     vertex.label.color="white",
     asp=FALSE,
     vertex.label.cex=0.4,
     edge.width=edgeweights,
     edge.arrow.mode=0,
     main="Co-occurrence modules", vertex.color=V(mst)$color)




#######################################
####Analisis all distance matrices#####
#######################################
#For this distance matrices between peptides in each module should be available

Compute_heatmap = function(Distance, Name, Subset_info, correlations){
        Out =   "" #Output path
        breaksList = seq(0, 1, by = 0.1)
        read.table(Distance, row.names=1, skip=1) -> Distance
        lapply(rownames(Distance), function(x){ str_split(x,":")[[1]][1] } ) %>% unlist() -> Names
        rownames(Distance) = Names
        colnames(Distance) = rownames(Distance)
        png(file= paste(Out, Name, "_similarity.png" ) )
        pheatmap::pheatmap(Distance, color=viridis::viridis(10), breaks = breaksList, main = Name)
        dev.off()
        
        Distance2 = Distance[upper.tri(Distance)]
        Similarity_score = mean(Distance2)
        Similarity_deviation = sd(Distance2)
        
        Data_LLD %>% select(colnames(Distance)) %>% cor() -> Corr_matrix
        test = hclust(dist(Corr_matrix))
        ORDER = rownames(Corr_matrix)[test$order]
        as.data.frame(Corr_matrix)[ORDER,ORDER] -> Corr_matrix ; Distance[ORDER, ORDER] -> Distance

        New_matrix = 1- Distance
        New_matrix[lower.tri(New_matrix)] <- Corr_matrix[lower.tri(Corr_matrix)]
        #New_matrix %>% drop_na() -> New_matrix
        pdf(file= paste(Out, Name, "_simVScorr.pdf" ) )
        corrplot::corrplot(as.matrix(New_matrix), diag=FALSE, tl.col="black")
        dev.off()
        
        #Get P-value
        vegan::mantel( as.dist(Corr_matrix) , as.dist(1-Distance), method = "spearman", permutations = 2000, na.rm = TRUE) -> Correlation_result
        r = Correlation_result$statistic
        p = Correlation_result$signif
        Result = tibble(Mean_similarity = Similarity_score, Sd_similarity=Similarity_deviation, r_mantel=r, P_mantel=p)
        return(Result)
        
}

files <- list.files(path="~/Desktop/Clusters", pattern="Distance_*", full.names=TRUE, recursive=FALSE)
read_tsv("Cluster_AtLeast10.tsv") -> Groups #Cluster belonging

Similarity_cluster = tibble()
for (Fi in files){
        basename(tools::file_path_sans_ext(Fi) ) -> N
        Groups %>% filter(Cluster == as.numeric(str_split(N, "_")[[1]][2]) ) ->Subset_info
        Compute_heatmap(Fi, N, Subset_info, correlations) -> D
        print(paste(N, D))
        rbind(Similarity_cluster, mutate(D, Cluster = str_split(N, "_")[[1]][2]) ) -> Similarity_cluster
}

Groups %>% filter(Cluster == 1) %>% select(Probe,  Taxa, Description) %>% print(n=99)
Similarity_cluster %>% filter(P_mantel>0.05)
write_tsv(Similarity_cluster,path = "Similarity_scores.tsv")




