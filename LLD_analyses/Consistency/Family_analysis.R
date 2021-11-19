library(tidyverse)
library(viridis)
library(ggforce)

Covariates_file =""
RDS_matrix = ""
Covariates_family =""

Prepare = F
if (Prepare == T){
        Cov = read_tsv(Covariates_file) #ID , sample_id
        Data = readRDS(RDS_matrix) 

        Data %>% filter(grepl("32_", ID)) -> Data
        Data$ID = str_remove(Data$ID, '"')
        Data$ID = as.vector(sapply(Data$ID, FUN= function(x){ str_split(x, "_")[[1]][2] } ))
        Data$sample_id = Data$ID
        select(Data, -ID) -> Data
        print(Data)
        
        Cov %>% filter(sample_id %in% Data$sample_id) -> Cov
        left_join(Cov, Data, by = "sample_id") -> DD

        DD %>% filter( grepl("LL", ID)) -> DD
        
		
        GoNL_info = read_tsv(Covariates_family
        colnames(GoNL_info)[1] = "ID"
        left_join(GoNL_info, DD) -> DD 
        print(DD)
        saveRDS(DD, file = "Data_analysis_family.rds")
} else{
	print("Reading")
        DD = readRDS(file ="Data_analysis_family.rds")
}


remove_l = function(x){
	n = c()
	for (m in x){
		z = substr(m,1,(nchar(m)-1))
		n = c(n,z)
	}
	return(n)
}
print("Preparing family info")
Family = remove_l(DD$WGS_id)
DD %>% mutate(Fam = Family) -> DD


Names = c("WGS_id", "relatives GoNL", "sample_id",  "Sex",  "Age", "plate_id", "Disease_status", "Fam")

#Pamily's distances

f <- function (i, j, dist_obj) {
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (i >= 1) & (j >= 1) & (i > j) & (i <= n) & (j <= n)
  k <- (2 * n - j) * (j - 1) / 2 + (i - j)
  k[!valid] <- NA_real_
  k
  }

Prevalence_f= function(x, Thresh = 0.1){
	V = sum(x, na.rm= T)/length(x)
	return( V >= Thresh)
}
print("Computing antibody prevalence")
AB_profile = select(DD, -c("ID", Names)) -> AB_profile
#Keep_ab = apply(AB_profile, 2, FUN = Prevalence_f) #Appplied on samples with two time points
#Keep_ab = colnames(AB_profile)[Keep_ab]

#AB_profile %>% select(Keep_ab) -> AB_profile
print("Computing distance matrix")
vegan::vegdist(AB_profile, method="jaccard", na.rm = T) -> Distance_matrix
print(DD %>% select(WGS_id, Fam))
DD %>% group_by(Fam) %>% summarise(N = n()) %>% group_by(N) %>% summarise(N_f = n() )  %>% print()


Distance_extraction =  function(DD, DM=Distance_matrix){
	DD %>% group_by(Fam) %>% summarise(N = n()) %>% filter(N == 3) -> Keep_fam
	Distances_common_tibble = tibble()
	#print("Substracting pairwise distances")
	No_trio = c()
	for (People in Keep_fam$Fam){
		#Take each person from the families that we want to check
		Look_up = which(DD$Fam == People) #Get all participants form the same family
		if (! length(Look_up) == 3){ No_trio = c(No_trio, People) ; next }
		k = f(Look_up[1], Look_up[2], DM) #Get comparison 1/2 in distance matrix
		k2 = f(Look_up[1], Look_up[3], DM)
		k3 = f(Look_up[2], Look_up[3], DM)
	
		ID1 = substr( DD$WGS_id[Look_up[1]],  nchar(DD$WGS_id[Look_up[1]]), nchar(DD$WGS_id[Look_up[1]]) )
		ID2 = substr(DD$WGS_id[Look_up[2]], nchar(DD$WGS_id[Look_up[2]]), nchar(DD$WGS_id[Look_up[2]]) )
		ID3 = substr(DD$WGS_id[Look_up[3]], nchar(DD$WGS_id[Look_up[3]]), nchar(DD$WGS_id[Look_up[3]]) )

		if (is.na(k))  { k  = f(Look_up[2], Look_up[1], DM) }
		if (is.na(k2)) { k2 = f(Look_up[3], Look_up[1], DM) } 
		if (is.na(k3)) { k3 = f(Look_up[3], Look_up[2], DM) } 	
	
		Distance1 = Distance_matrix[k]
		Distance2 = Distance_matrix[k2]
		Distance3 = Distance_matrix[k3]
		id1 = paste(c(ID1, ID2), collapse = "-") ; id2 = paste(c(ID1, ID3), collapse = "-") ; id3 = paste(c(ID2, ID3), collapse = "-")	
	
		New_t = tibble(ID = c(People, People, People), Pair= c("1-vs-2", "1-vs-3", "2-vs-3") , k = c(k,k2,k3), Distance = c(Distance1, Distance2, Distance3), relation = c(id1,id2, id3 ))
		Distances_common_tibble = rbind(Distances_common_tibble, New_t)
	}

	Distances_common_tibble %>% drop_na() -> Distances_common_tibble
	Distances_common_tibble  %>% group_by(Pair) %>% summarise(M = mean(Distance)) -> Stat
	return(list(Distances_common_tibble, Stat) )
}

Permutations = 2000
Distances =  Distance_extraction(DD)
Distances_common_tibble = Distances[[1]]
Stats = Distances[[2]]

set.seed(99)
Permuted_stats = tibble()
for (perm in seq(Permutations)) {
	#Permute family labels
	DD_p = DD
	DD_p$Fam = DD_p$Fam[ sample(nrow(DD_p))]
	DE = Distance_extraction(DD_p)
	rbind(Permuted_stats, DE[[2]]) -> Permuted_stats
}
print(Stats)
print(Permuted_stats)

print("Comparison to permutation H0")
perm_P = tibble()
for (C in c("1-vs-2", "1-vs-3", "2-vs-3") ){
	Stats %>% filter(Pair == C) %>% select(M) %>% as_vector() %>% as.vector() -> f
	Permuted_stats %>% filter(Pair == C) %>% select(M) %>% as_vector() %>% as.vector() -> distribution_0	
	
	length(distribution_0[distribution_0<=f])/Permutations -> P
	rbind(perm_P, tibble(Comparison = C, `P-value`=P)) -> perm_P
}
print(perm_P)




Distances_common_tibble %>% group_by(ID) %>% summarise(n())
Distances_common = Distances_common_tibble$k
Common = Distance_matrix[Distances_common]
Not_common =  Distance_matrix[-Distances_common]
print(dim(Distance_matrix))


wilcox.test(Common, Not_common) -> Stat_w
P_diff = Stat_w$p.value
print( paste(c("Wilcox distance related-unrealated: ", as.character(P_diff)), collapse="") )

Distance_tibble = rbind( tibble(Distance = Common, Family = T), tibble(Distance = Not_common, Family = F) )
#Distance_tibble %>% mutate(Family = ifelse(Family=="mum-dad", "Mother to Father", ifelse(Family=="dad-child", "Father to Offspring", ifelse(Family == "mum-child", "Mother to Offspring", Family)))) %>% mutate(Family = factor(Family, levels=c("unrelated", "Mother to Father", "Father to Offspring", "Mother to Offspring" ))) -> plot_input


Distance_tibble  %>% ggplot(aes(x=Family, y = Distance)) + geom_violin(aes(fill=Family)) + theme_bw() +scale_color_manual("Family relationship" ,values = c("#440154FF", "#3B528BFF", "#21908CFF") ) + ylab("Distance (Jaccard)")  -> Distance_box
Distance_tibble %>% ggplot(aes(x=Distance, fill = Family)) + geom_density(alpha=0.5) +theme_bw() + scale_color_manual("Family relationship" ,values = c("#440154FF", "#3B528BFF", "#21908CFF") ) -> Distance_distribution




Dad_child = Distance_matrix[ filter(Distances_common_tibble, relation == "a-c" )$k ]
Mum_child = Distance_matrix[ filter(Distances_common_tibble, relation == "b-c" )$k ]
mum_dad = Distance_matrix[ filter(Distances_common_tibble, relation == "a-b" )$k ]

Distance_tibble = rbind(tibble(Distance = Dad_child, Family = "Father to Offspring"), tibble(Distance = Mum_child, Family = "Mother to Offspring"), tibble(Distance = mum_dad, Family = "Mother to Father"), tibble(Distance = Not_common, Family = "unrelated") )

Distance_tibble %>% mutate(Family = factor(Family, levels = c("unrelated", "Mother to Father", "Father to Offspring", "Mother to Offspring"))) -> Distance_tibble

Distance_tibble  %>% ggplot(aes(x=Family, y = Distance)) + geom_boxplot(outlier.shape = NA,aes(col=Family)) + geom_sina(alpha = 0.2) + theme_bw() + scale_color_manual("Family relationship", values = c("#440154FF", "#35B779FF", "#31688EFF", "#80CDC1") ) + ylab("Distance (Jaccard)") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank()) -> Distance_box

#Distance_tibble %>% ggplot(aes(x=Family, y = Distance, col = Family)) + geom_boxplot() +geom_point() +theme_bw() -> Distance_box
lm( Distance ~ Family, Distance_tibble) -> M1
TUKEY <- TukeyHSD(x=aov(M1), 'Family', conf.level=0.95)
as.data.frame(TUKEY$Family) %>% rownames_to_column("Comparison") -> TUKEY





