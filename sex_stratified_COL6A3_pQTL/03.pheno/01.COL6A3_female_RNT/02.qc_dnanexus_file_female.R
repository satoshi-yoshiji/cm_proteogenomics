library(tidyverse)
library(data.table)
library(magrittr)
library(readxl)
library(RNOmni)

setwd('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_female_RNT/')

# pheno <- fread('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3/02.dnanexus/COL6A3_pQTL_pheno_participant.tsv', fill = F)
# colnames(pheno)
# nrow(pheno) # 488127
# 
# # Create a named vector for the Genetic Principal Components columns renaming
# PC_names <- setNames(paste0("p22009_a", 1:20), paste0("PC", 1:20))
# 
# # Use the select function to rename and keep only the specified columns
# pheno <- pheno %>%
#   select(
#     FID = eid,
#     IID = eid,
#     Sex = p31,
#     Genetic_sex = p22001,
#     Age_at_recruitment = p21022,
#     UK_Biobank_assessment_centre = p54_i0,
#     Time_blood_sample_collected = p3166_i0_a0,
#     Sex_chromosome_aneuploidy = p22019,
#     Genotype_measurement_batch = p22000,
#     Genetic_ethnic_grouping = p22006,
#     !!!PC_names
#   )
# 
# # transform Sex and Genetic_sex columns to binary (0/1)
# pheno %<>% mutate(Sex = ifelse(Sex == 'Female', 0, 1)) %>% mutate(Genetic_sex = ifelse(Genetic_sex == 'Female', 0, 1)) # Female = 0 , Male = 1
# 
# # one hot encoding for UK_Biobank_assessment_centre
# # Step 1: One-hot encoding
# pheno_encoded <- pheno %>%
#   mutate(Flag = 1) %>%
#   spread(key = UK_Biobank_assessment_centre, value = Flag, fill = 0)
# 
# # Step 2: Rename columns
# new_colnames <- c("centre" %>% paste0(seq_len(length(unique(pheno$UK_Biobank_assessment_centre)))))
# names(pheno_encoded)[which(names(pheno_encoded) %in% unique(pheno$UK_Biobank_assessment_centre))] <- new_colnames
# 
# # Genotype_measurement_batch
# #  transform the Genotype_measurement_batch column such that"Batch_" = 0 and "UKBiLEVEAX_" = 1,
# pheno_encoded %<>%
#   mutate(Genotype_measurement_batch = case_when(
#     startsWith(Genotype_measurement_batch, "Batch_") ~ 0,
#     startsWith(Genotype_measurement_batch, "UKBiLEVEAX_") ~ 1,
#     TRUE ~ NA_integer_
#   ))
# 
# # QC pheno
# pheno_qc <- pheno_encoded %>% 
#   filter(Sex == Genetic_sex) %>% select(-Genetic_sex) %>%  # sex mismatch
#   filter(Sex_chromosome_aneuploidy != 'Yes') %>% select(-Sex_chromosome_aneuploidy) %>% # sex chromosome aneuploidy
#   filter(Genetic_ethnic_grouping == 'Caucasian') %>% select(-Genetic_ethnic_grouping) # Caucasian only
#   
# # check  
# head(pheno_qc)  
# 
# ###################  
# # merge protein level data (outcome)
# ###################  
# protlev <- fread('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3/01.olink_pheno.tsv', fill = T)
# protlev %<>% transmute(IID = eid, Protein_level = result, Batch, Processing_StartDate)
# protlev %<>% mutate(Batch = as.numeric(Batch))
# protlev %<>% filter(Batch >=1) %>% filter(Batch <= 6) # exclude Batch 0 (pilot) and Batch 7 (COVID)
# 
# nrow(protlev) # 44831
# 
# # one-hot encoding for Batch
# protlev_encoded <- protlev %>% 
#   mutate(Batch_olink1 = case_when(Batch == 1 ~ 1, TRUE ~ 0),
#          Batch_olink2 = case_when(Batch == 2 ~ 1, TRUE ~ 0),
#          Batch_olink3 = case_when(Batch == 3 ~ 1, TRUE ~ 0),
#          Batch_olink4 = case_when(Batch == 4 ~ 1, TRUE ~ 0),
#          Batch_olink5 = case_when(Batch == 5 ~ 1, TRUE ~ 0),
#          Batch_olink6 = case_when(Batch == 6 ~ 1, TRUE ~ 0)#,
#          #Batch7 = case_when(Batch == 7 ~ 1, TRUE ~ 0)
#          ) 
# protlev_encoded %<>% select(-Batch) 
# 
# # combine 
# combine <- inner_join(protlev_encoded, pheno_qc, by = 'IID')
# combine %<>% relocate(FID, IID, Sex, Age_at_recruitment, Time_blood_sample_collected)
# combine %<>% drop_na() # 36872
# 
# # add "time_between_collect_process"
# combine %<>% mutate(
#   Processing_StartDate = ymd(Processing_StartDate),
#   Time_blood_sample_collected = ymd_hms(Time_blood_sample_collected),
#   Time_between_collect_process = as.numeric(difftime(Processing_StartDate, Time_blood_sample_collected, units = "days")))
# combine %<>% select(-Processing_StartDate, -Time_blood_sample_collected) 
# 
# # check 
# nrow(combine) # 36872
# 
# # finally, create all required covariates
# # For the discovery cohort, association models included the following covariates: age, age2, sex,
# # age*sex, age2*sex, batch, UKB centre, UKB genetic array, time between blood sampling and
# # measurement and the first 20 genetic principal components (PCs).
# combine %<>% mutate(Age2 = Age_at_recruitment*Age_at_recruitment, 
#                     AgeSex = Age_at_recruitment*Sex, 
#                     Age2Sex = Age_at_recruitment*Age_at_recruitment*Sex)
# 
# 
# # relocate columns
# combine %<>% relocate(FID, IID, Protein_level, 
#                       Age_at_recruitment,
#                       Age2,
#                       Sex, 
#                       AgeSex,
#                       Age2Sex,
#                       Time_between_collect_process)
# 
# # remmove duplicated individuals
# combine_unique <- combine %>% filter(!duplicated(IID))
# nrow(combine_unique) # 36670
# 
# #check final colnames
# colnames(combine_unique)


# read 02.pheno_covariate.tsv
combine_unique <- fread('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_RNT/02.pheno_covariate.tsv') # Note this is not normalized yet

# Sex = 0 (female) only
combine_unique_female <- combine_unique %>% filter(Sex == 0)
nrow(combine_unique_female) # 19770

# RNT
combine_unique_female$Protein_level <- RNOmni::RankNorm(combine_unique_female$Protein_level)

# save
write_tsv(combine_unique_female, file = '02.pheno_covariate_female.tsv')
write_tsv(combine_unique_female %>% select(FID, IID, Protein_level) , file = '02.pheno_female.tsv')
write_tsv(combine_unique_female %>% select(-Protein_level) , file = '02.covariate_female.tsv')

# save individuals to be kept for step 1 and 2 analyses
keep_female0 <- combine_unique_female %>% select(FID, IID) 
exclude <- fread('/scratch/richards/satoshi.yoshiji/08.UKB/excluded_individuals/w27449_2023-04-25.tsv')  
#define Not in operator function
`%!in%` <- Negate(`%in%`)
keep_female_excluded <- keep_female0 %>% filter(!IID %in% exclude$V1) 

# write
write_tsv(keep_female_excluded, file = '02.keep_female.tsv', col_names = F)



