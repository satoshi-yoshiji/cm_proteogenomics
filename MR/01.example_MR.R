library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)
library(MRPRESSO)

wd <- "/scratch/yoshiji/18.obesity_cardiometab/00.BMI_no_MHC/"
setwd(wd)

output_dir <- paste0(wd, "output/")
system(paste0('mkdir -p ', output_dir, 'harmonized/'))
system(paste0('mkdir -p ', output_dir, 'or/'))
system(paste0('mkdir -p ', output_dir, 'hetero/'))
system(paste0('mkdir -p ', output_dir, 'pleio/'))
system(paste0('mkdir -p ', output_dir, 'steiger/'))
system(paste0('mkdir -p ', output_dir, 'mrpresso/'))
system(paste0('mkdir -p ', output_dir, 'mrpresso_add/'))


# read the outcome GWAS
outcome_path <- '/home/yoshiji/scratch/11.pQTL/03.pQTL_decode2021/11196_31_COL6A3_Collagen_alpha_3_VI_.txt.gz'
protname <- '11196_31_COL6A3_Collagen_alpha_3_VI'

#args <- commandArgs(trailingOnly=TRUE) #input UNZIPPED file
#outcome_path <- args[1]
#protname <- args[2]

outcome_GWAS <- vroom(outcome_path)
outcome_GWAS %<>% dplyr::rename(SNP = rsids) #IEUGWAS format

#exposure
exp_path <- '/home/yoshiji/scratch/09.proMR/01.exposure/BMI/BMI_noMHC_clumped.tsv'
exp_dat <- read_exposure_data(filename = exp_path,
                              sep='\t',
                              snp_col='SNP',
                              beta_col = 'beta.exposure',
                              se_col = 'se.exposure',
                              effect_allele_col = 'effect_allele.exposure',
                              other_allele_col = 'other_allele.exposure',
                              eaf_col = 'eaf.exposure',
                              pval_col = 'pval.exposure')

# # To find proxies, you may use either LDlinkR, snappy, or ieugwasr

#######################
# common 
# downstream procedures are the same from MR with or without the use of proxies
#######################

# read the outcome
outcome_GWAS %<>% dplyr::filter(SNP %in% exp_dat$SNP) 
formatted_outcome <- format_data(outcome_GWAS, snps = exp_dat$SNP, type="outcome", snp_col="SNP", beta_col="Beta", se_col = "SE", eaf_col = "ImpMAF", 
                             effect_allele_col = "effectAllele", other_allele_col = "otherAllele", pval_col = "Pval", chr_col = "Chrom", pos_col = "Pos")

# harmonize
exp_dat_outcome <-harmonise_data(exposure_dat=exp_dat, outcome_dat=formatted_outcome)
# exclude rare variants
exp_dat_outcome %<>% filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999))
exp_data_outcome_name <- paste0(output_dir,"harmonized/", protname, ".harmonized.txt")
write.table(exp_dat_outcome, file=exp_data_outcome_name, sep = '\t', quote = F, row.names = F)

# mr result
mr_results <- mr(exp_dat_outcome)
#mr_name <- paste0(output_dir, "mr/", protname, ".mr.txt")
#write.table(mr_results, file=mr_name, sep = '\t', quote = F)

# odds ratio
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write.table(OR, file=OR_name, sep = '\t', quote = F, row.names = F)

# # scatter plot
# pdf_name <- paste0(output_dir, "pdf/", protname, ".scatter.pdf")
# pdf(pdf_name)
# mr_scatter_plot(mr_results, exp_dat_outcome)[[1]] + xlim(0,1)
# dev.off()

# horizontal pleiotropy
tryCatch({ 
pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
pleio_name <- paste0(output_dir,"pleio/", protname, ".pleio.txt")
write.table(pleio_res, file=pleio_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# hetero test
tryCatch({ 
hetero_res <- mr_heterogeneity(exp_dat_outcome)
hetero_res$isquared <- abs(100*(hetero_res$Q - hetero_res$Q_df)/hetero_res$Q)  # I2 = 100%×(Q - df)/Q
hetero_name <- paste0(output_dir,"hetero/", protname, ".hetero.txt")
write.table(hetero_res, file=hetero_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# steiger
tryCatch({ 
exp_dat_outcome$samplesize.exposure <- 681275
exp_dat_outcome$samplesize.outcome <- 35559  
steiger <- directionality_test(exp_dat_outcome)
steiger_name <- paste0(output_dir, "steiger/", protname, ".steiger.txt")
write.table(steiger, file=steiger_name, sep = '\t', quote = F, row.names = F)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# MR-PRESSO
tryCatch({ 
  mrpresso_res <- mr_presso(BetaExposure = "beta.exposure",
                            BetaOutcome = "beta.outcome",
                            SdOutcome = "se.outcome",
                            SdExposure = "se.exposure",
                            OUTLIERtest = TRUE,
                            DISTORTIONtest = TRUE,
                            data = exp_dat_outcome,
                            NbDistribution = 15000,
                            SignifThreshold = 0.05)
  
  mrpresso_name <- paste0(output_dir, "mrpresso/", protname, ".mrpresso.txt")
  
  # add global_rss, global_pval, distortion_indices, distortion_coef, distortion_pval
  mrpresso_df <- as.data.frame(mrpresso_res$`Main MR results`)
  mrpresso_df %<>% mutate(#global_rss = mrpresso_res$`MR-PRESSO results`$`Global Test`$RSSobs,
    global_pval =  mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue,
    #distorition_indices = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,
    #distortion_coef = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`,
    distortion_pval = mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue)
  
  write_tsv(mrpresso_df, file =  mrpresso_name)
  
  # additional information
  mrpresso_global_rss <- mrpresso_res$`MR-PRESSO results`$`Global Test`$RSSobs
  mrpresso_global_pval <- mrpresso_res$`MR-PRESSO results`$`Global Test`$Pvalue
  mrpresso_distortion_indices <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  mrpresso_distortion_coef <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`
  mrpresso_distortion_pval <-  mrpresso_res$`MR-PRESSO results`$`Distortion Test`$Pvalue
  
  mrpresso_add_res <- data.frame(mrpresso_global_rss = mrpresso_global_rss,
                                 mrpresso_global_pval = mrpresso_global_pval,
                                 mrpresso_distortion_indices = mrpresso_distortion_indices,
                                 mrpresso_distortion_coef = mrpresso_distortion_coef,
                                 mrpresso_distortion_pval = mrpresso_distortion_pval
  )
  
  mrpresso_add_name <- paste0(output_dir, "mrpresso_add/", protname, ".mrpresso_add.txt")
  write.table(mrpresso_add_res, file=mrpresso_add_name, sep = '\t', quote = F, row.names = F)
  
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
