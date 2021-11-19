library(tidyverse)

readRDS("") -> df #RDS of presence/absence 
read_tsv("Data/Probes_to_keep.txt") -> filter #List of pepetides that need to stay



sapply(colnames(df), function(x){ gsub('"', '',x) } ) -> df_names
colnames(df) = df_names

print(dim(df))
df %>% select(c("ID", filter$value)) -> df
print(dim(df))

#Write New matrix that will be used for all other analyses
write_rds(df, "Data/Immuno_matrix_postSelection.rds")
write_tsv(df, "Data/Immuno_matrix_postSelection.tsv")


