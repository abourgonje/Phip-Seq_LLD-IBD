library(tidyverse)

Merge = function(ID1, ID2){
	Path_files = "../Heritability/2.GREML/Phenotypes/"
	File1 = paste(c(Path_files, ID1), collapse="")
	File2 = paste(c(Path_files, ID2), collapse="")	

	Pheno1 = read_tsv(File1, col_names= F)
	Pheno2 = read_tsv(File2, col_names= F)
	
	colnames(Pheno1)[3] = "Pheno1"
	colnames(Pheno2)[3] = "Pheno2"
	left_join(Pheno1, Pheno2) -> Combined

	N = paste(c(ID1, ID2), collapse = "-")
	Name = paste(c("Phenotypes/", N), collapse="")

	write_tsv(Combined, Name, col_names= F)
	return(Name)
}



List_probes = read_tsv("Probes_to_run.txt", col_names=F)
Not_do = read_tsv("Probes_already_run.txt", col_names=F)
List_probes %>% filter(! X1 %in% Not_do$X1) -> List_probes
script = "Parser.py"
N = length(List_probes$X1)
Counter = 0
read_tsv("GeneticCorrelation_LLD_20heritabilityIncomplete.tsv") -> Incomplete
Phenos_Available = read_tsv("Phenos_available", col_names = F)
print(Phenos_Available)

Done = c(Incomplete$ID)
for (P in seq(N)){
	Probe = List_probes$X1[P]
	for (P2 in seq(P+1:N)){
		Probe2 = List_probes$X1[P2]
		if (Probe == Probe2){ next }
		paste(sort(c(Probe, Probe2)), collapse="") -> ID
		#if(ID %in% Done){ next }
		if (ID %in% Phenos_Available$X1){ next }
		Counter = Counter + 1
		Name = Merge(Probe, Probe2)
		command = paste(c("python", script, Name),collapse=" ")
		#system(command)	
}
}


print(Counter)


