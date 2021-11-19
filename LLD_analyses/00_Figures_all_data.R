#Produce figures and stats using all the dataset before filtering by prevalence/sequence identity
library(tidyverse)

Uncompressed_data_matrix = "" #Presence/absence matrix to be read in tsv format
Compressed_data_matrix = "" #Same as uncompressed_data_matrix but saved in an RDS file
Covariates = "" #Covariate file

#1 read
Do = F #Reading takes a loooong time if not compressed. Adjust this if you have Compressed_data_matrix
if (Do == T){ 
Data = read_tsv(Uncompressed_data_matrix)
saveRDS(Data, file = Compressed_data_matrix)
}else { Data = readRDS(file = Compressed_data_matrix ) }
#The prefix indicates the dataset (IBD/CeD/LLD), we remove CeD cohort
Data %>% filter(! grepl("34_", ID) ) -> Data
print("All data, dimensions")
print(dim(Data))
print("Dimensions LLD")
Data %>% filter( grepl("32_", ID) ) %>% dim() %>% print()
#2. Count number of antibodies per participant and per cohort
Count_AB_N = function(Data){
        dplyr::select(Data, -ID) %>% apply(1, function(x) { sum(x,na.rm = T) } )  %>% as_vector() -> Number_ab
        sapply(Data$ID, function(x){ str_split(x, "_")[[1]][1] } )  %>% as_vector() -> Cohort 
        tibble(ID = Data$ID, Cohort = Cohort, Probe_n = Number_ab) -> N_AB_tibble
        return(N_AB_tibble)
}

Count_AB_N(Data) -> N_AB_tibble
print("Summary of antibody number")
summary(N_AB_tibble$Probe_n)



#Plot distribution of number of probes
mutate(N_AB_tibble, Cohort = ifelse(Cohort == "32", "LLD", "1000IBD" )) %>% ggplot(aes(x=log10(Probe_n), fill = Cohort ) )+ geom_density(alpha =0.5) + theme_bw() -> Fig 
ggsave(filename = "Results/Exploration_allProbes/Distribution-number.pdf" , plot = Fig)
N_AB_tibble %>% filter(Cohort == "32") %>% ggplot(aes(x=log10(Probe_n)) )+ geom_density() + theme_bw() -> Fig2
ggsave(filename = "Results/Exploration_allProbes/Distribution-number_LLD.pdf" , plot = Fig)


print("Samples with really low number of probes")
N_AB_tibble %>% arrange(Probe_n) %>% print()
# A tibble: 2,281 x 3
#   ID            Cohort Probe_n
#   <chr>         <chr>    <dbl>
# 1 32_3260408512 32           3
# 2 32_3260414603 32          12
# 3 32_SF00086694 32          28
# 4 32_SF00374694 32          75
# 5 32_SF00170602 32         142
# 6 33_SF02185603 33         179
# 7 32_SF03619693 32         196
# 8 33_82499320   33         205
# 9 33_1042345307 33         211
#10 32_SF00115588 32         213




#match with covariates
sapply(N_AB_tibble$ID, function(x){ str_split(x, "_")[[1]][2] } )  %>% as_vector() -> Names
N_AB_tibble$ID = Names
read_tsv(Covariate)   %>%  filter(sample_id %in% N_AB_tibble$ID ) %>% mutate(Sample_name = ID ,ID = sample_id) %>% select(- sample_id) %>% drop_na() -> Covs

left_join(N_AB_tibble, Covs, by="ID") -> N_AB_tibble
N_AB_tibble %>% group_by(Cohort) %>% summarise( Mean = mean(Probe_n), SD = sd(Probe_n), Max = max(Probe_n), Min = min(Probe_n) ) %>% print()
print("Do each cohort have a different number of probes")
lm( log10(Probe_n) ~ Cohort + Age + Sex, N_AB_tibble ) -> M_comparison_nprobe
summary(M_comparison_nprobe) %>% print() # IBD has on average fewer +probes. Age has a positive effect on number of probes. Sex has not effect



#2.Count prevalence per antibody
Prevalence = function(Data){
        dplyr::select(Data, -ID) %>% apply(2, function(x){ sum(x, na.rm = T) }  )  %>% as.data.frame() %>% rownames_to_column("Probe") -> Data2
        colnames(Data2) = c("Probe", "Number of samples")
        Data2 %>% mutate(`Percentage of samples` = `Number of samples`/dim(Data)[1]) -> Data2
	Data2 %>% group_by(`Number of samples`) %>% summarise(N = n() ) %>% mutate(N+1) -> Data3
        return(list(Data2, Data3))
}
Prevalence_all = Prevalence(Data)
Prevalence_LLD = Prevalence(filter(Data, grepl("32", ID) ) ) 
Prevalence_IBD = Prevalence(filter(Data, grepl("33", ID) ) ) 

Prevalence_all[[1]] %>% mutate(Cohort = "Both") -> P_all ; Prevalence_all[[2]] %>% mutate(Cohort = "Both") -> N_all
Prevalence_IBD[[1]] %>% mutate(Cohort = "IBD") -> P_IBD	; Prevalence_IBD[[2]] %>% mutate(Cohort = "IBD") -> N_IBD
Prevalence_LLD[[1]] %>% mutate(Cohort = "LLD") -> P_LLD	; Prevalence_LLD[[2]] %>% mutate(Cohort = "LLD") -> N_LLD
rbind(rbind(P_all, P_IBD), P_LLD) -> Prevalence



#3. Barplot 

Prevalence -> Input_figure1B

ggplot(Prevalence, aes(x=`Number of samples`)) + geom_histogram(aes(fill= Cohort), alpha=0.6 ) + theme_bw() + scale_y_log10() -> Number_plot
ggsave(filename = "Results/Exploration_allProbes/Number_shared.pdf" , plot = Number_plot)
Prevalence %>% filter(Cohort=="LLD") %>% ggplot(aes(x=`Number of samples`)) + geom_histogram() + theme_bw() + scale_y_log10() -> Number_plot_LLD
ggsave(filename = "Results/Exploration_allProbes/Number_shared_LLD.pdf" , plot = Number_plot_LLD)



for (C in unique(Prevalence$Cohort)){
	print(paste(c("Probes with prevalence of at least 1", C), collapse=" "))
	Prevalence %>% filter(Cohort == C) %>%  group_by(`Number of samples` >= 1) %>% summarise( N = n() ) %>% print()
	print(paste(c("Quantiles of number of people per probe", C), collapse=" ")) 
	quantile(filter(Prevalence, Cohort==C)$`Number of samples`, probs = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1 )) %>% print()
	#filter(Prevalence, Cohort==C) %>% filter(`Number of samples` > 1000) %>% arrange(`Number of samples`) %>% print()
	ecdf(filter(Prevalence, Cohort==C)$`Number of samples`) -> P
	P(1000) - 1  %>% print()
}

