# Step 1 MR

# load libraries
library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)
library(MRPRESSO)

# Define working directory and set output directories
wd <- "path-to-your-wd"
setwd(wd)
output_dir <- paste0(wd, "output/")
dirs <- c("harmonized", "or", "hetero", "pleio", "steiger", "mrpresso", "mrpresso_add")
lapply(dirs, function(dir) system(paste0('mkdir -p ', output_dir, dir, '/')))

args <- commandArgs(trailingOnly=TRUE) 
outcome_path <- args[1]
protname <- args[2]

# Read the outcome GWAS
outcome_GWAS <- vroom(outcome_path) %>% dplyr::rename(SNP = rsids) 

# Read the exposure data
exp_path <- 'path-to-your-exposure-data' # BMI_noMHC_clumped.tsv
exp_dat <- read_exposure_data(
  filename = exp_path,
  sep = '\t',
  snp_col = 'SNP',
  beta_col = 'beta.exposure',
  se_col = 'se.exposure',
  effect_allele_col = 'effect_allele.exposure',
  other_allele_col = 'other_allele.exposure',
  eaf_col = 'eaf.exposure',
  pval_col = 'pval.exposure'
)

# Filter outcome data to match exposure SNPs
outcome_GWAS <- outcome_GWAS %>% dplyr::filter(SNP %in% exp_dat$SNP)
formatted_outcome <- format_data(outcome_GWAS, snps = exp_dat$SNP, type = "outcome", 
                                 snp_col = "SNP", beta_col = "Beta", se_col = "SE", eaf_col = "ImpMAF", 
                                 effect_allele_col = "effectAllele", other_allele_col = "otherAllele", 
                                 pval_col = "Pval", chr_col = "Chrom", pos_col = "Pos")

# For proxy search, you may use snappy v1.0 (https://gitlab.com/richards-lab/vince.forgetta/snappy/-/blob/master/snappy)

# Harmonize data
exp_dat_outcome <- harmonise_data(exposure_dat = exp_dat, outcome_dat = formatted_outcome) %>%
  filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999))

# Write harmonized data
exp_data_outcome_name <- paste0(output_dir, "harmonized/", protname, ".harmonized.txt")
write.table(exp_dat_outcome, file = exp_data_outcome_name, sep = '\t', quote = F, row.names = F)

# Run MR and write results
mr_results <- mr(exp_dat_outcome)
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write.table(OR, file = OR_name, sep = '\t', quote = F, row.names = F)

# Pleiotropy test
tryCatch({
  pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
  pleio_name <- paste0(output_dir, "pleio/", protname, ".pleio.txt")
  write.table(pleio_res, file = pleio_name, sep = '\t', quote = F, row.names = F)
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n") })

# Heterogeneity test
tryCatch({
  hetero_res <- mr_heterogeneity(exp_dat_outcome)
  hetero_res$isquared <- abs(100 * (hetero_res$Q - hetero_res$Q_df) / hetero_res$Q)
  hetero_name <- paste0(output_dir, "hetero/", protname, ".hetero.txt")
  write.table(hetero_res, file = hetero_name, sep = '\t', quote = F, row.names = F)
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n") })

# Steiger directionality test
tryCatch({
  exp_dat_outcome$samplesize.exposure <- 681275
  exp_dat_outcome$samplesize.outcome <- 35559
  steiger <- directionality_test(exp_dat_outcome)
  steiger_name <- paste0(output_dir, "steiger/", protname, ".steiger.txt")
  write.table(steiger, file = steiger_name, sep = '\t', quote = F, row.names = F)
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n") })

# MR-PRESSO test
tryCatch({
  mrpresso_res <- mr_presso(
    BetaExposure = "beta.exposure",
    BetaOutcome = "beta.outcome",
    SdOutcome = "se.outcome",
    SdExposure = "se.exposure",
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    data = exp_dat_outcome,
    NbDistribution = 15000,
    SignifThreshold = 0.05
  )
  saveRDS(mrpresso_res, file = paste0(output_dir, "mrpresso/", protname, ".mrpresso.RDS"))
  
  mrpresso_df <- as.data.frame(mrpresso_res$`Main MR results`) %>%
    mutate(global_pval = mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue,
           distortion_pval = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue)
  mrpresso_name <- paste0(output_dir, "mrpresso/", protname, ".mrpresso.txt")
  write_tsv(mrpresso_df, file = mrpresso_name)
  
  # Additional MR-PRESSO results
  mrpresso_add_res <- data.frame(
    mrpresso_global_rss = mrpresso_res$`MR-PRESSO results`$`Global Test`$RSSobs,
    mrpresso_global_pval = mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue,
    mrpresso_distortion_indices = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,
    mrpresso_distortion_coef = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`,
    mrpresso_distortion_pval = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue
  )
  mrpresso_add_name <- paste0(output_dir, "mrpresso_add/", protname, ".mrpresso_add.txt")
  write.table(mrpresso_add_res, file = mrpresso_add_name, sep = '\t', quote = F, row.names = F)
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n") })
