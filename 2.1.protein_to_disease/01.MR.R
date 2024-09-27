# Step 2 MR

# load libraries
library(tidyverse)
library(vroom)
library(magrittr)
library(TwoSampleMR)

output_dir <- paste0(wd, "output/")
dirs <- c("harmonized", "or", "hetero", "pleio", "steiger")
lapply(dirs, function(dir) system(paste0('mkdir -p ', output_dir, dir, '/')))

# Dynamic input from command line arguments
args <- commandArgs(trailingOnly = TRUE)
exp_path <- args[1] # path to cis-pQTL
protname <- args[2] # path to outcome

# Read exposure data
exp_dat <- read_exposure_data(
  filename = exp_path,
  sep = '\t',
  snp_col = 'variant',
  beta_col = 'beta_unadj',
  se_col = 'se',
  effect_allele_col = 'Amin',
  other_allele_col = 'Amaj',
  eaf_col = 'MAF',
  pval_col = 'pval'
)

# Use snappy and/or LDlinkR for proxy search

# Read and format outcome data
formatted_outcome <- format_data(outcome_GWAS, snps = exp_dat$SNP, type = "outcome",
                                 snp_col = "SNP", beta_col = "ES", se_col = "SE", eaf_col = "AF",
                                 effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "pval",
                                 chr_col = "seqnames", pos_col = "start", samplesize_col = "SS")
formatted_outcome$id.outcome <- 'outcome'

# Harmonize data
exp_dat_outcome <- harmonise_data(exposure_dat = exp_dat, outcome_dat = formatted_outcome)

# Write harmonized data
exp_data_outcome_name <- paste0(output_dir, "harmonized/", protname, ".harmonized.txt")
write_tsv(exp_dat_outcome, file = exp_data_outcome_name)

# MR results and odds ratio
mr_results <- mr(exp_dat_outcome)
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write_tsv(OR, file = OR_name)

# Pleiotropy test
pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
pleio_name <- paste0(output_dir, "pleio/", protname, ".pleio.txt")
write_tsv(pleio_res, file = pleio_name)

# Heterogeneity test with i-squared
tryCatch({
  res_single <- mr_singlesnp(exp_dat_outcome)
  res_Isq <- Isq(res_single$b, res_single$se)
  hetero_res <- mr_heterogeneity(exp_dat_outcome)
  hetero_res$i_squared <- res_Isq
  hetero_name <- paste0(output_dir, "hetero/", protname, ".hetero.txt")
  write_tsv(hetero_res, file = hetero_name)
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n") })

# Steiger test
exp_dat_outcome$samplesize.exposure <- 35559
steiger <- directionality_test(exp_dat_outcome)
steiger_name <- paste0(output_dir, "steiger/", protname, ".steiger.txt")
write_tsv(steiger, file = steiger_name)

# Final print message
print(paste0(protname, ": done"))
