# Load libraries
library(tidyverse)
library(vroom)
library(TwoSampleMR)
library(ieugwasr)

# Load BMI GWAS data
gwas <- vroom('path-to-wd/BMI/Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt')

# Remove SNPs in the MHC region (chr6:28,477,797-33,448,354, GRCh37)
mhc_region <- gwas %>%
  arrange(CHR, POS) %>%
  filter(CHR == 6, POS >= 28477797, POS <= 33448354)

mhcsnp <- mhc_region %>%
  select(SNP)

# Filter out MHC region SNPs from the full GWAS dataset
gwas_nomhc <- gwas %>%
  filter(!(SNP %in% mhcsnp$SNP))

# Format the GWAS data for exposure
exp_dat <- format_data(dat = gwas_nomhc,
                       type = 'exposure',
                       snp_col = 'SNP',
                       beta_col = 'BETA',
                       se_col = 'SE',
                       eaf_col = 'Freq_Tested_Allele_in_HRS',
                       effect_allele_col = 'Tested_Allele',
                       other_allele_col = 'Other_Allele',
                       pval_col = 'P',
                       samplesize_col = 'N')

# Perform clumping on the exposure data
exp_clump <- exp_dat %>%
  mutate(rsid = SNP, pval = pval.exposure) %>%
  ld_clump(plink_bin = "path-to-plink/plink", 
           bfile = "path-to-1KGEUR-ref", 
           clump_p = 5e-8)

# Save the clumped results to a file
write.table(exp_clump, file = 'BMI_noMHC_clumped.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
