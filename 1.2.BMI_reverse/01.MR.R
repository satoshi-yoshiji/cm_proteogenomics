# Reverse MR for step 1

# Load necessary libraries
library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)

# Set working directory
wd <- "path-to-your-wd"
setwd(wd)

# Create output directories
output_dir <- paste0(wd, "output/")
dirs <- c("harmonized", "or", "hetero", "pleio", "steiger")
lapply(dirs, function(dir) system(paste0('mkdir -p ', output_dir, dir, '/')))

# Read outcome GWAS
outcome_path <- 'path-to-ukb-b-19953.tsv' # Use BMI GWAS of UKB for more variants, including cis-pQTL of proteins
outcome_GWAS <- vroom(outcome_path) %>%
  dplyr::rename(SNP = ID)  # Rename for IEUGWAS format

# Dynamic input for exposure data via command line arguments
args <- commandArgs(trailingOnly = TRUE)
exp_path <- args[1]
protname <- args[2]

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

# Format outcome data for harmonization
formatted_outcome <- format_data(
  outcome_GWAS, snps = exp_dat$SNP, type = "outcome",
  snp_col = "SNP", beta_col = "ES", se_col = "SE", eaf_col = "AF",
  effect_allele_col = "ALT", other_allele_col = "REF",
  pval_col = "pval", chr_col = "seqnames", pos_col = "start"
)
formatted_outcome$id.outcome <- 'outcome'

# For proxy search, you may use snappy v1.0 (https://gitlab.com/richards-lab/vince.forgetta/snappy/-/blob/master/snappy)

# Harmonize exposure and outcome data
exp_dat_outcome <- harmonise_data(exposure_dat = exp_dat, outcome_dat = formatted_outcome) %>%
  filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999))

# Save harmonized data
exp_data_outcome_name <- paste0(output_dir, "harmonized/", protname, ".harmonized.txt")
write_tsv(exp_dat_outcome, file = exp_data_outcome_name)

# Run MR analysis and save results
mr_results <- mr(exp_dat_outcome)

# Generate and save odds ratios
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write_tsv(OR, file = OR_name)

# Conduct horizontal pleiotropy test and save results
pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
pleio_name <- paste0(output_dir, "pleio/", protname, ".pleio.txt")
write_tsv(pleio_res, file = pleio_name)

# Conduct heterogeneity test and calculate I-squared
tryCatch({
  res_single <- mr_singlesnp(exp_dat_outcome)
  res_Isq <- Isq(res_single$b, res_single$se)
  hetero_res <- mr_heterogeneity(exp_dat_outcome)
  hetero_res$i_squared <- res_Isq
  hetero_name <- paste0(output_dir, "hetero/", protname, ".hetero.txt")
  write_tsv(hetero_res, file = hetero_name)
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
})

# Conduct Steiger directionality test and save results
# Ensure that the sample size variables are defined beforehand
exp_dat_outcome$samplesize.exposure <- samplesize.exposure  # Define exposure sample size
exp_dat_outcome$samplesize.outcome <- samplesize.outcome    # Define outcome sample size
steiger <- directionality_test(exp_dat_outcome)
steiger_name <- paste0(output_dir, "steiger/", protname, ".steiger.txt")
write_tsv(steiger, file = steiger_name)

# Print completion message
print(paste0(protname, ": done"))
