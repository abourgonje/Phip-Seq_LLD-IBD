#Compute distance between samples with two timepoints

args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(viridis)
if (length(args) == 1){ 
	if (args[[1]] == "Repeat"){
		Prepare = T
	} else { Prepare = F }
} else { Prepare = F}

Covariates_file = ""
RDS_matrix = ""

if (Prepare == T){
        Cov = read_tsv(Covariates_file) #ID , sample_id
        Data = readRDS(RDS_matrix) 
	
        Data %>% filter(grepl("32_", ID)) -> Data
        Data$ID = str_remove(Data$ID, '"')
        Data$ID = as.vector(sapply(Data$ID, FUN= function(x){ str_split(x, "_")[[1]][2] } ))
        Data$sample_id = Data$ID
        select(Data, -ID) -> Data
        
        Cov %>% filter(sample_id %in% Data$sample_id) -> Cov
        left_join(Cov, Data, by = "sample_id") -> DD

        DD %>% filter(grepl("LL", ID)) -> DD
        DD %>% mutate(Time_point = ifelse(grepl("_F", ID), "FA", "BL")) -> DD
        DD$ID =  as.vector(sapply(DD$ID, FUN= function(x){ str_split(x, "_F")[[1]][1] } ))
        DD_FA = filter(DD, Time_point == "FA")
        DD_BL = filter(DD, Time_point == "BL")


        DD_FA %>% filter(ID %in% DD_BL$ID) -> DD_FA
        DD %>% filter( ID %in% DD_FA$ID) -> DD

        saveRDS(DD, file = "Data_analysis_FA.rds")
}else{
        # Restore the object
        DD = readRDS(file = "Data_analysis_FA.rds")
}



DD %>% group_by(Time_point) %>% summarise(n())

Names = c("sample_id", "Sex", "Age", "plate_id", "Time_point","Disease_status")

#Participant's distances

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

AB_profile = select(DD, -c("ID", Names)) -> AB_profile
#Keep_ab = apply(AB_profile, 2, FUN = Prevalence_f) #Appplied on samples with two time points
#Keep_ab = colnames(AB_profile)[Keep_ab]

#AB_profile %>% select(Keep_ab) -> AB_profile

Get_distances = function(DD, Distance_matrix){
	Distances_common_tibble = tibble()
	for (People in unique(DD$ID)){
        	Look_up = which(DD$ID == People)
	        k = f(Look_up[1], Look_up[2], Distance_matrix)
        	if (is.na(k)) { k = f(Look_up[2], Look_up[1], Distance_matrix) }
	        if (is.na(k)) { next }
        	Distance = Distance_matrix[k]
	        Distances_common_tibble = rbind(Distances_common_tibble, tibble(ID = People, k= k, Distance= Distance_matrix[k] ))
	}

	Distances_common = Distances_common_tibble$k
	Common = Distance_matrix[Distances_common]
	Not_common =  Distance_matrix[-Distances_common]
	return( list(Common, Not_common, mean(Common)  ) )
}
vegan::vegdist(AB_profile, method="jaccard", na.rm = T) -> Distance_matrix
Permutations = 2000
Get_distances(DD, Distance_matrix) -> Not_permuted
Common = Not_permuted[[1]] ;  Not_common = Not_permuted[[2]]
t = Not_permuted[[3]]
Null_distribution = c()
for (p in Permutations){
	DD_p = DD
	DD_p$ID = DD_p$ID[sample(nrow(DD_p)) ]
	Get_distances(DD_p, Distance_matrix) -> permuted
	Null_distribution = c(Null_distribution, permuted[[3]] )	
}




wilcox.test(Common, Not_common) -> Stat_w
P_diff = Stat_w$p.value
print( paste(c("Wilcox distance related-unrealated: ", as.character(P_diff)), collapse="") )


p = length(Null_distribution[Null_distribution<=t])/Permutations
print( paste(c("Pvalue to null distribution of means from permutations: ", as.character(p)), collapse="") )



Distance_tibble = rbind( tibble(Distance = Common, Longitudinal = T), tibble(Distance = Not_common, Longitudinal = F) )

saveRDS(Distance_tibble,"/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Immuno_markers/Input_figures/Fig1D.rds")
q()
Distance_tibble %>% ggplot(aes(x=Longitudinal, y = Distance)) + geom_violin(aes(fill=Longitudinal )) +theme_bw() + scale_fill_viridis_d("Longitudinal samples") + ylab("Distance (Jaccard)") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) -> Distance_box
Distance_tibble %>% ggplot(aes(x=Distance, fill = Longitudinal)) + geom_density(alpha=0.5) +theme_bw() + scale_colour_viridis_d("viridis") -> Distance_distribution



#What factors determine change?
Distances_common_tibble %>% arrange(desc(Distance))
left_join(Distances_common_tibble, select(DD, c(ID, Sex, Age) )  , by="ID") %>% drop_na() -> Distances_common_tibble
Distances_common_tibble %>% mutate(Age_n = Age-min(Age) ) -> Distances_common_tibble
print(Distances_common_tibble)
lm(Distance ~ Sex+ Age_n, Distances_common_tibble) -> M
summary(M)




###Analysis per antibody

AB_profile = select(DD, -c("ID", Names)) -> AB_profile
Keep_ab = apply(AB_profile, 2, FUN = Prevalence_f, Thresh=0.05)
Keep_ab = colnames(AB_profile)[Keep_ab]




#Do_consistency = T


if (Prepare  == T){
        #By antibody consistency
        Antibody_table = tibble()
        for ( Antibody in Keep_ab  ){
                if (Antibody %in% c(Names, "ID")){ next }
                DD %>% select(c(Antibody)) %>% as.vector() %>% as_vector() -> AB

                DD %>% select(Names) %>% mutate(Anb = AB) -> Antibody_data
        
                Antibody_data1 =   filter(Antibody_data, Time_point == "BL") %>% arrange(sample_id)
                Antibody_data2 = filter(Antibody_data, Time_point == "FA") %>% arrange(sample_id)

                Antibody_data1 %>% mutate( BL = Antibody_data1$Anb , FA = Antibody_data2$Anb  ) %>% select(-c("Time_point", Anb)) ->Antibody_data_wide
                #Antibody_data %>% spread("Time_point", Antibody) -> Antibody_data_wide
        
        
                NA_rate_BL = sum(is.na(Antibody_data_wide$BL)) / dim(Antibody_data_wide)[1]
                NA_rate_FA = sum(is.na(Antibody_data_wide$FA)) / dim(Antibody_data_wide)[1]

                Antibody_data_wide %>% drop_na() -> Antibody_data_wide
        
                Prevalence_BL = sum(Antibody_data_wide$BL , na.rm= T) / dim(Antibody_data_wide)[1]
                Prevalence_FA = sum( Antibody_data_wide$FA , na.rm= T) / dim(Antibody_data_wide)[1]

                #C = cor(Antibody_data_wide$BL, Antibody_data_wide$FA, method = "pearson" ) #"p"
                Consistency = sum( Antibody_data_wide$BL == Antibody_data_wide$FA ) / dim(Antibody_data_wide)[1] 
		Change_mean = mean( Antibody_data_wide$FA - Antibody_data_wide$BL )

                N_T = tibble("Antibody" = Antibody, "Consistency" = Consistency, "Mean_change" = Change_mean, "Prevalence_BL" = Prevalence_BL, "Prevalence_FA" = Prevalence_FA, "NA_rate_BL" =NA_rate_BL, "NA_rate_FA"=NA_rate_FA )
		        
                Antibody_table = rbind(Antibody_table, N_T)     
}
        saveRDS(Antibody_table, file = "Antibody_consistency.rds")

} else { Antibody_table =  readRDS(file ="Antibody_consistency.rds")  }


#read_tsv("Prevalence_antibodies_total.tsv") -> Prevalence_info
#colnames(Prevalence_info) = c("Antibody","Prevalence_LLD")
#read_csv("../Data/df_info_AT_v2 2.csv") -> Meta_data #order, full name
#Meta_data %>% mutate(Antibody = order) %>% select(Antibody, `full name`) -> Meta_data
#left_join(Antibody_table, Meta_data, by="Antibody") -> Antibody_table
#left_join(Antibody_table, Prevalence_info, by="Antibody") -> Antibody_table
#Antibody_table %>% arrange(Consistency) %>% print()
Antibody_table %>% ggplot(aes(x=Consistency )) + geom_density() + theme_bw() -> Plot_consistency
Antibody_table %>% group_by(Mean_change >= 0 ) %>% summarise(n()) %>% print()




