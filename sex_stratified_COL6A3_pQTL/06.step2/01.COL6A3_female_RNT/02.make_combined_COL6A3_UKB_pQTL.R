library(vroom)
library(tidyverse)
library(magrittr)

setwd('/home/richards/satoshi.yoshiji/scratch/30.pQTL_regenie/06.step2/01.COL6A3_female_RNT')

per_chr_pqtl_list <- list.files('output', full.names = T)
# Filtering the list to retain only those that end with ".regenie"
filtered_list <- grep("\\.regenie$", per_chr_pqtl_list, value = TRUE)

# Extract chromosome numbers
chromosome_numbers <- as.integer(gsub(".*chr([0-9]+)_.*", "\\1", filtered_list))

# Order the list by chromosome numbers
ordered_list <- filtered_list[order(chromosome_numbers)]

# Print the ordered list
ordered_list

combined_gwas <- data.frame()
for(i in ordered_list){
  ith_pqtl <- vroom(i)
  combined_gwas <- rbind(combined_gwas, ith_pqtl)
}

# save
write_tsv(combined_gwas, file = 'COL6A3.combined.Protein_level.regenie.gz')
