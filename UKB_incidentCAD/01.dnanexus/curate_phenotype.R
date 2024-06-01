library(tidyverse)
library(magrittr)
library(vroom)

setwd("/scratch/richards/satoshi.yoshiji/18.obesity_cardiometab/42.UKB_incidentCAD_obs")

pheno <- vroom('phenotype.tsv')

# remove [] and '' from ICD10 column (p41270)
pheno$p41270 <- gsub("\\[|\\]|'", "", pheno$p41270)

# count items in the ICD10 diagnosis column
pheno$p41270_item_count <- sapply(strsplit(pheno$p41270, ",\\s*"), length)
pheno %<>% relocate(p41270_item_count, .after = 'p41270')

pheno1 <- pheno %>% mutate(p41280_item_count = rowSums(!is.na(select(., starts_with("p41280_")))))
pheno1 %<>% relocate(p41280_item_count, .after = 'p41270_item_count')
# Display the head of the updated data frame to verify the new column
head(pheno1)


pheno1 <- pheno[3,]
