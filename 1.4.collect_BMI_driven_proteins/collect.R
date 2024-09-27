# Load necessary libraries
library(vroom)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggrepel)

# Set working directory
setwd("path-to-wd")

# Import the BMI_noMHC-to-protein result with significance, heterogeneity, and pleiotropy tests
nsum <- vroom('path-to-results-of-1.1.BMI_to_proteins', delim = '\t')

# Get causal estimates - IVW
nsum2 <- nsum %>%
  filter(method == 'Inverse variance weighted') %>%
  group_by(protein) %>%
  arrange(pval) %>%
  filter(!duplicated(protein)) %>%
  ungroup() %>%
  mutate(shortname = sapply(strsplit(protein, "_"), function(x) x[3]),
         seqid = paste0(sapply(strsplit(protein, "_"), function(x) x[1]),
                        '_',
                        sapply(strsplit(protein, "_"), function(x) x[2]))) %>%
  select(protein, shortname, seqid, method, nsnp, b, se, pval, lo_ci, up_ci)

# Weighted median, mode, and Egger results
nsum2_median <- nsum %>%
  filter(method == 'Weighted median') %>%
  group_by(protein) %>%
  arrange(pval) %>%
  filter(!duplicated(protein)) %>%
  ungroup() %>%
  select(protein, b, se, pval, lo_ci, up_ci) %>%
  rename(b.median = b, se.median = se, lo_ci.median = lo_ci, up_ci.median = up_ci, pval.median = pval)

nsum2_mode <- nsum %>%
  filter(method == 'Weighted mode') %>%
  group_by(protein) %>%
  arrange(pval) %>%
  filter(!duplicated(protein)) %>%
  ungroup() %>%
  select(protein, b, se, pval, lo_ci, up_ci) %>%
  rename(b.mode = b, se.mode = se, lo_ci.mode = lo_ci, up_ci.mode = up_ci, pval.mode = pval)

nsum2_egger <- nsum %>%
  filter(method == 'MR Egger') %>%
  group_by(protein) %>%
  arrange(pval) %>%
  filter(!duplicated(protein)) %>%
  ungroup() %>%
  select(protein, b, se, pval, lo_ci, up_ci) %>%
  rename(b.egger = b, se.egger = se, lo_ci.egger = lo_ci, up_ci.egger = up_ci, pval.egger = pval)

# Get heterogeneity and pleiotropy results for IVW
hetero <- vroom('path-to-results-of-1.1.BMI_to_protein-heterogeneity-results.txt', delim = '\t') %>%
  filter(!is.na(Q_pval)) %>%
  select(protein, Q, Q_df, Q_pval, isquared)

pleio <- vroom('/path-to-results-of-1.1.BMI_to_protein-pleiotropy-results.txt', delim = '\t') %>%
  filter(!is.na(pval)) %>%
  select(protein, egger_intercept, se, pval) %>%
  rename(se.egger_intercept = se, pval.egger_intercept = pval)

# Join all results (IVW, median, mode, Egger, hetero, pleio)
njoin <- nsum2 %>%
  left_join(hetero, by = 'protein') %>%
  left_join(pleio, by = 'protein') %>%
  left_join(nsum2_median, by = 'protein') %>%
  left_join(nsum2_mode, by = 'protein') %>%
  left_join(nsum2_egger, by = 'protein') %>%
  group_by(protein) %>%
  arrange(pval) %>%
  filter(!duplicated(protein)) %>%
  ungroup() %>%
  separate(protein, into = c('seqid1', 'seqid2', 'shortname'), sep = '_', remove = FALSE) %>%
  mutate(seqid = paste(seqid1, seqid2, sep = '_')) %>%
  relocate(seqid, .after = 'protein') %>%
  select(-seqid1, -seqid2) %>%
  mutate(protein = paste(shortname, seqid, sep = '.'))

# Annotate significance, heterogeneity, and pleiotropy results
njoin <- njoin %>%
  mutate(Significance = factor(case_when(pval < 1e-5 ~ 'Pass', TRUE ~ 'No'), levels = c('No', 'Pass')),
         Heterogeneity_test = case_when(isquared < 50 ~ 'Pass', TRUE ~ 'Fail'),
         MREgger_intercept_test = case_when(pval.egger_intercept > 0.05 ~ 'Pass', TRUE ~ 'Fail'))

# Incorporate MR-PRESSO results
mrpresso_res_collect <- read_tsv('path-to-results-of-1.1.BMI_to_protein-mrpresso-results.tsv')
njoin <- left_join(njoin, mrpresso_res_collect, by = 'seqid')

# Outlier robust estimate test with IVW, weighted median, Egger, and MR-PRESSO
njoin <- njoin %>%
  mutate(Outlier_robust_estimate_test = case_when(
    (b > 0 & b.median > 0 & b.egger > 0 & is.na(mrpresso.causal_estimate.corrected)) |
      (b < 0 & b.median < 0 & b.egger < 0 & is.na(mrpresso.causal_estimate.corrected)) |
      (b > 0 & b.median > 0 & b.egger > 0 & !is.na(mrpresso.causal_estimate.corrected) & mrpresso.causal_estimate.corrected > 0) |
      (b < 0 & b.median < 0 & b.egger < 0 & !is.na(mrpresso.causal_estimate.corrected) & mrpresso.causal_estimate.corrected < 0) ~ 'Pass',
    TRUE ~ 'Fail'))

# Incorporate reverse causation and additional tests
reverse <- vroom('path-to-results-of-1.2.BMI_reverse', delim = '\t')
reverse_join <- reverse %>%
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) %>%
  group_by(protein) %>%
  arrange(pval) %>%
  filter(!duplicated(protein)) %>%
  ungroup()

reverse_join <- reverse_join %>%
  left_join(hetero, by = 'protein') %>%
  left_join(pleio, by = 'protein') %>%
  mutate(Significance = factor(case_when(pval < 0.05/nrow(reverse_join) ~ 'Pass', TRUE ~ 'No')),
         Heterogeneity_test = case_when(isquared < 0.5 ~ 'Pass', TRUE ~ 'Fail'),
         MREgger_intercept_test = case_when(pval.egger_intercept >= 0.5 ~ 'Pass', TRUE ~ 'Fail')) %>%
  mutate(Sensitivity_test = case_when(
    (Heterogeneity_test == 'Pass' | is.na(Heterogeneity_test)) &
      (MREgger_intercept_test == 'Pass' | is.na(MREgger_intercept_test)) ~ 'Pass',
    TRUE ~ 'Fail'))

# Final filter and output
njoin_clean <- njoin %>%
  select(protein, seqid, shortname, method, nsnp, b, se, pval, lo_ci, up_ci, Q, Q_df, Q_pval, isquared,
         b.median, se.median, pval.median, lo_ci.median, up_ci.median, b.egger, se.egger, pval.egger, lo_ci.egger, up_ci.egger,
         mrpresso.global_pval, mrpresso.causal_estimate.corrected, bodyfatpercentage.b, bodyfatpercentage.pval,
         Significance, Heterogeneity_test, MREgger_intercept_test, Outlier_robust_estimate_test)

write_tsv(njoin_clean, 'aptamer_summary.tsv')