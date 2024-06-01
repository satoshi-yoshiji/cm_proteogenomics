library(tidyverse)
library(vroom)
library(magrittr)
library(readxl)
#library(RNOmni) # RNT will be conducted in "02.qc_dnanexus_file.R"

setwd("/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_RNT/")

coding <-vroom('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/data/coding143.tsv')
olink_assay <- vroom('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/data/olink_assay.dat')
olink_assay_version <- vroom('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/data/olink_assay_version.dat')
olink_batch_number <- vroom('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/data/olink_batch_number.dat')
olink_data <- vroom('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/data/olink_data.txt')
olink_processing_start_date <- vroom('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/data/olink_processing_start_date.dat')
plate_participant <- vroom('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/data/UKBPPP-plate_participant.csv', delim = ',')

# data
olink_data %<>% select(eid, protein_id, result)
# coding
coding %<>% rename(protein_id = coding)

################
# restrict olink_data to COL6A3 
# protein_id == 637
# meaning = COL6A3;Collagen alpha-3(VI) chain
################
olink_data %<>% filter(protein_id == '637')
#olink_data$result <- RankNorm(olink_data$result)
#coding %<>% filter(meaning == 'COL6A3;Collagen alpha-3(VI) chain')
#olink_assay %<>% filter(Assay == 'COL6A3')

################
# curation starts here
################

# olink_df
olink_df <- left_join(olink_data, coding, by = 'protein_id')
olink_df %<>% separate(meaning, into = c('symbol', 'full'), sep = ';', remove = F)

# olink_assay (Assay   UniProt Panel)
olink_df2 <- left_join(olink_df, olink_assay %>% rename(symbol = Assay), by = 'symbol')

# assay version  (Panel_Lot_Nr Assay  Assay_Version)
olink_df3 <- left_join(olink_df2, olink_assay_version %>% rename(symbol = 'Assay'), by = 'symbol')

# plate_participant (eid p30901_i0)
plate_participant %<>% dplyr::select(eid, p30901_i0) %>% rename(PlateID = p30901_i0) %>% filter(!is.na(PlateID))
olink_df4 <- left_join(olink_df3, plate_participant, by = 'eid')

# olink_batch_number (PlateID Batch)
olink_df5 <- left_join(olink_df4, olink_batch_number, by = 'PlateID')

# olink_processing_start_date (PlateID       Panel           Processing_StartDate)
olink_df6 <- left_join(olink_df5, olink_processing_start_date, by = c('PlateID', 'Panel'))
# length(unique(olink_df6$eid)) [1] 52749
# length(unique(olink_df6$protein_id)) [1] 1463
# length(unique(olink_df6$symbol)) [1] 1463
# length(unique(olink_df6$UniProt)) [1] 1463
# unique(olink_df6$Panel) #"Cardiometabolic" "Neurology"  "Inflammation" "Oncology"  
# unique(olink_df6$Panel_Lot_Nr)) # "B04414" "B04413" "B04412" "B04411"
# unique(olink_df6$Assay_Version) all are one
# length(unique(olink_df6$PlateID)) 633 -> Not used
# unique(olink_df6$Batch) 4  7  6  5  3  2  1 NA  0;  length(unique(olink_df6$Batch)) = 9

# save
write_tsv(olink_df6 %>% dplyr::select(-meaning, -Panel_Lot_Nr, -Assay_Version, -PlateID), file = '01.olink_pheno.tsv')


