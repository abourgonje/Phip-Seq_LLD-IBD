library(tidyverse)
library(Cairo)
#MarkerName Allele1 Allele2 Weight Zscore P-value Direction Probe


#1. Prepare Significant/or nearly significant associations

#../GWAS/agilent_242450.assoc.logistic
#files <- list.files(path="../GWAS", pattern="*.logistic", full.names=TRUE, recursive=FALSE)
#Total = tibble()
#for (f in files){
#	Base_name = str_split(f, "/")[[1]][3]
#	N = str_split(Base_name, "\\.")[[1]][1]
#	File = read.table(f) %>% as_tibble() %>% mutate(Phenotype = N) -> File
#	rbind(Total, File) -> Total
#}
#read_tsv("Check", col_types= c("c","n", "n", "n", "c", "c") ) -> C
#print(C)

Do = T
if (Do == T){
	Total =  read.table("Most_significant_GWAS.txt", header=T) %>% as_tibble()
	#Trash = read.table("Trash_GWAS.txt", header=T) %>% as_tibble()
	#rbind(Total, Trash) -> Total
	lapply(Total$MarkerName, function(x){ str_split(x, ":")[[1]]  } ) -> Splitted
	Chr = sapply(Splitted, function(x){ x[1]} )
	BP = sapply(Splitted, function(x){ x[2]} )
	Total %>% mutate(CHR = as.factor(Chr), BP=as.numeric(BP), P=`P.value`) ->  Total
	saveRDS(Total, file = "GWAS.rds")
} else { 
	Total = readRDS("GWAS.rds")
	Total %>% mutate(CHR = as.factor(CHR)) ->  Total
}

#print(Total)

# CHR          SNP         BP
#For GWAS
Study_wide = 5.67e-11
Genome_Wide = 5e-8

Total %>% filter(P < Genome_Wide) %>% select(-P) -> Total2
write_tsv(Total2, "Highly_significant.txt")

Total %>% filter(! grepl("\\?",Direction)) -> Total

#Get numbers of chormosomic positions
#Order factors
Order = as.character(c(seq(22)))

Total$CHR <- factor(Total$CHR, levels = Order)

Chr_length = tibble(CHR= as.factor(seq(1:22)), L=c(248956422, 242193529, 198295559,190214555, 181538259,
                                                                          170805979, 159345973 ,145138636, 138394717, 133797422, 135086622,
                                                                          133275309, 114364328, 107043718, 101991189, 90338345,
                                                                          83257441, 80373285, 58617616, 64444167,46709983 ,50818468) )
L2 = c(0, cumsum(Chr_length$L))
Chr_length %>% mutate(L2=L2[1:22]) -> Chr_length
left_join(Total, Chr_length, by="CHR") %>% mutate(BPcum = BP + L2) -> Total






#Total$BPcum <- NA
#s <- 0
#nbp <- c()

#For each chromosome get the maximum of the position, and add it up to the total
#for (i in unique(sort(as.factor(Total$CHR)))){
#  nbp[i] <- max(Total[Total$CHR == i,]$BP)
#  Cum =tibble( Total[Total$CHR == i,"BP"] + s )
#  colnames(Cum) = "BP"
#  Total %>% mutate(BPcum = ifelse(CHR == i, Cum$BP, BPcum)) -> Total
  #Total[Total$CHR == i,"BPcum"] <- Total[Total$CHR == i,"BP"] + s
  #print(Total[Total$CHR == i,"BPcum"])
#  s <- s + nbp[i]
#}
print(Total)

#Prepare axis
axis.set <- Total %>%
  group_by(CHR) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(Total$P)))) + 3
nCHR = length(unique(Total$CHR))

THRESH = 1 #1e-5


BAR = tibble(CHR =  axis.set$CHR,BP= axis.set$center, P= -log10(THRESH) )
#BAR = tibble(whole.genome$BPcum, chr =  whole.genome$chr, P= -log10(THRESH))
#Relevelinag
BAR$CHR <- factor(BAR$CHR, levels = Order)
axis.set$CHR = factor(axis.set$CHR, levels = Order)
 

library(ggrepel)
manhplot =  ggplot() +
  geom_point(data = Total, aes(x = BPcum, y = -log10(P), color = CHR, size = -log10(P)), alpha = 0.5) +
  geom_hline(yintercept = -log10(Study_wide), color = "#dbc0d9", size=1.5) +
  geom_hline(yintercept = -log10(Genome_Wide), color = "#9fa4cd", linetype = "dashed", size=1.5) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(-log10(4), ylim)) +
  scale_color_manual(values = rep(c("#a5a5a5", "#8c8c8d"), nCHR)) +
  scale_fill_manual(values = rep(c("#a5a5a5", "#8c8c8d"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL,
       y = "-log10(p)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  ) #+ ylim(0, NA)

#tiff( "Manh.tiff", height = 13.6, width = 13.6, units = "in", res=800)
#dev.off()
ggsave(plot = manhplot, filename = "Manh.tiff", height = 13.6, width = 13.6, units = "in")






