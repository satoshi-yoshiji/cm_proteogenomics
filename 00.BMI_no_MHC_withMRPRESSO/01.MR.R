library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)
library(MRPRESSO)
#library(data.table) # for copy() function
#library(LDlinkR)

wd <- "/scratch/yoshiji/18.obesity_cardiometab/00.BMI_no_MHC/" #!!! don't forget the slash (/) at the end of the full path
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
outcome_path <- '/scratch/richards/public/decode_proteomics_2021/4337_49_CRP_CRP.txt.gz'
protname <- '4337_49_CRP_CRP'
# outcome_path <- '/home/yoshiji/scratch/11.pQTL/03.pQTL_decode2021/11196_31_COL6A3_Collagen_alpha_3_VI_.txt.gz'
# protname <- '11196_31_COL6A3_Collagen_alpha_3_VI'

#args <- commandArgs(trailingOnly=TRUE) #input UNZIPPED file
#outcome_path <- args[1]
#protname <- args[2]
outcome_GWAS <- vroom(outcome_path)
outcome_GWAS %<>% dplyr::rename(SNP = rsids) #IEUGWAS format


# exp_dat <- read_exposure_data(filename = exp_path,
#                               sep='\t',
#                               snp_col='',
#                               beta_col = 'beta_unadj',
#                               se_col = 'se',
#                               effect_allele_col = 'Amin',
#                               other_allele_col = 'Amaj',
#                               eaf_col = 'MAF',
#                               pval_col = 'pval')

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

# # find proxy
# if(file.exists(paste0('/scratch/richards/satoshi.yoshiji/18.obesity_cardiometab/00.fatpercentage_to_deCODE_proxy/snappy/', protname, '/snappyout.txt'))){
# proxy_snp_cov <- vroom(paste0('/scratch/richards/satoshi.yoshiji/18.obesity_cardiometab/00.fatpercentage_to_deCODE_proxy/snappy/', protname, '/snappyout.txt'), col_names = F)
# colnames(proxy_snp_cov)<-c("SNP", "SNP_out", "r2", "V4", "match_check")
# 
# if(length(proxy_snp_cov) == 0){
# proxy_snp_cov_filtered_snp <-proxy_snp_cov # do nothing if no proxy SNP
# } else{
# proxy_snp_cov_filtered_snp <-proxy_snp_cov %>% filter(SNP %in% exp_dat$SNP) # exposure = bodyfat; outcome = protein
# # table(proxy_snp_cov_filtered_snp$match_check)
# }
# 
# #######################
# # do proxy matching only if there are snappy results
# #######################
# if('proxy' %in% proxy_snp_cov_filtered_snp$match_check ){ 
#   
# MR_LeadSNPs_proxy <- left_join(x = exp_dat, y = proxy_snp_cov_filtered_snp[,c(1,2)], by = 'SNP')
# MR_LeadSNPs_proxy_1 <-copy(MR_LeadSNPs_proxy)
# MR_LeadSNPs_proxy_1$SNP_out <- ifelse(MR_LeadSNPs_proxy_1$SNP_out == ".", MR_LeadSNPs_proxy_1$SNP, MR_LeadSNPs_proxy_1$SNP_out)
# 
# ## merge the lead snps from exposure outcome GWAS by SNP_out
# MR_LeadSNPs_proxy_2<-merge(MR_LeadSNPs_proxy_1, outcome_GWAS, by.x="SNP_out", by.y = "SNP", all.x=TRUE)
# 
# # since there is not allele frequencies from meta-analysis, use allele1 frequency from european using dbsnp
# MR_LeadSNPs_proxy_2_proxy_set<-MR_LeadSNPs_proxy_2%>%filter(SNP_out != SNP)
# 
# allele_matching_table<-data.frame(
#   orgi_snp=character(),
#   prox_snp=character(),
#   orgi_Allele=character(),
#   cor_proxy_allele=character())
# 
# re <- "\\(([^()]+)\\)"  # get a character inside the bracket
# for (i in 1:nrow(MR_LeadSNPs_proxy_2_proxy_set)){
#   tryCatch({
#   original_snp=MR_LeadSNPs_proxy_2_proxy_set$SNP[i]
#   proxy_snp=MR_LeadSNPs_proxy_2_proxy_set$SNP_out[i]
#     ##get original and proxy snp information
#     ld_pair_result<-LDpair(original_snp, proxy_snp, token="57241ca07ea6", api_root = 'https://ldlink.nci.nih.gov/LDlinkRest2/')
#     ld_pair_result <- ld_pair_result %>% filter(r2 >= 0.8) # somethimes rare SNPs can cause discordant results between snappy and LDpair. Here I ensured r2 >= 0.8 (e.g., s558161750 rs576003299)
#     #extract corresponding allele pair 1
#     allele_match1<-unlist(str_split(str_split(ld_pair_result$corr_alleles,", ")[[1]][1], "-"))
#     orginal_snp1=allele_match1[grepl(original_snp,allele_match1)] # grep original_snp (rsid) in allele_match1 (e.g., "rs35350651("   ")"             "rs10774625(A)")
#     cor_proxy_snp1=allele_match1[grepl(proxy_snp,allele_match1)]
#     orginal_allele1<-gsub(re, "\\1", str_extract_all(orginal_snp1, re)[[1]])  # get a character inside the bracket
#     cor_proxy_allele1<-gsub(re, "\\1", str_extract_all(cor_proxy_snp1, re)[[1]])
#     #extract corresponding allele pair 1 (SY note :2?)
#     allele_match2<-unlist(str_split(str_split(ld_pair_result$corr_alleles,", ")[[1]][2], "-"))
#     orginal_snp2=allele_match2[grepl(original_snp,allele_match2)]
#     cor_proxy_snp2=allele_match2[grepl(proxy_snp,allele_match2)]
#     orginal_allele2<-gsub(re, "\\1", str_extract_all(orginal_snp2, re)[[1]])
#     cor_proxy_allele2<-gsub(re, "\\1", str_extract_all(cor_proxy_snp2, re)[[1]])
#     #create a temp data_frame to save this line of result
#     temp_dat<-data.frame(orgi_snp=original_snp, prox_snp=proxy_snp, orgi_Allele=c(orginal_allele1,orginal_allele2), cor_proxy_allele=c(cor_proxy_allele1, cor_proxy_allele2))
#     allele_matching_table<-rbind(allele_matching_table, temp_dat)
#     
#     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#     }
# 
# MR_LeadSNPs_proxy_2_set1_matching<-MR_LeadSNPs_proxy_2%>%filter(SNP_out == SNP)##the set no need to be changed (matching snps)
# MR_LeadSNPs_proxy_2_setNA<-MR_LeadSNPs_proxy_2%>%filter(is.na(SNP_out) ==T)##the set no need to be changed (no matching snps) 
# MR_LeadSNPs_proxy_2_proxy_set<-MR_LeadSNPs_proxy_2%>%filter(is.na(SNP_out) ==F & SNP_out != SNP)## set need to be changed (with proxy snps)
# 
# ##fix A1,A2 (no allele will be matched for the snps that have different allele from MR and 1000genome - multiallelic/indel)
# MR_LeadSNPs_proxy_2_proxy_set1<-MR_LeadSNPs_proxy_2_proxy_set%>%rowwise()%>%
#   mutate(effect_allele.exposure = ifelse(effect_allele.exposure %in% allele_matching_table$orgi_Allele[allele_matching_table$orgi_snp == SNP], allele_matching_table$cor_proxy_allele[allele_matching_table$orgi_snp == SNP & allele_matching_table$orgi_Allele== effect_allele.exposure], NA),
#          other_allele.exposure = ifelse(other_allele.exposure %in% allele_matching_table$orgi_Allele[allele_matching_table$orgi_snp == SNP], allele_matching_table$cor_proxy_allele[allele_matching_table$orgi_snp == SNP & allele_matching_table$orgi_Allele== other_allele.exposure], NA))
# MR_LeadSNPs_proxy_2_proxy_set2<-MR_LeadSNPs_proxy_2_proxy_set1%>%filter(is.na(effect_allele.exposure)!=T & is.na(other_allele.exposure)!=T)
# 
# MR_LeadSNPs_proxy_3<-rbind(MR_LeadSNPs_proxy_2_set1_matching, MR_LeadSNPs_proxy_2_proxy_set2,MR_LeadSNPs_proxy_2_setNA)
# 
# ## Redefine variables as to reflect increasing exposure levels, which are associated with lower metabolites: 
# # Flipping Alleles and allele frequency
# MR_LeadSNPs_proxy_4<-copy(MR_LeadSNPs_proxy_3)
# MR_LeadSNPs_proxy_4$eaf.exposure <- ifelse(MR_LeadSNPs_proxy_3$beta.exposure > 0, MR_LeadSNPs_proxy_3$eaf.exposure, 1-MR_LeadSNPs_proxy_3$eaf.exposure)
# MR_LeadSNPs_proxy_4$effect_allele.exposure <- ifelse(MR_LeadSNPs_proxy_3$beta.exposure > 0, as.character(MR_LeadSNPs_proxy_3$effect_allele.exposure), as.character(MR_LeadSNPs_proxy_3$other_allele.exposure))
# MR_LeadSNPs_proxy_4$other_allele.exposure <- ifelse(MR_LeadSNPs_proxy_3$beta.exposure > 0, as.character(MR_LeadSNPs_proxy_3$other_allele.exposure), as.character(MR_LeadSNPs_proxy_3$effect_allele.exposure))
# MR_LeadSNPs_proxy_4$beta.exposure <- ifelse(MR_LeadSNPs_proxy_3$beta.exposure > 0, MR_LeadSNPs_proxy_3$beta.exposure, (-1)*MR_LeadSNPs_proxy_3$beta.exposure)
# ## add a column for study name
# MR_LeadSNPs_proxy_4$Study<-rep("MR", nrow(MR_LeadSNPs_proxy_4))
# 
# # remove MHC region
# MR_LeadSNPs_proxy_5 <- copy(MR_LeadSNPs_proxy_4)
# MR_LeadSNPs_proxy_5 %<>% filter(!(Chrom == 'chr6' & Pos >= 28477797 & Pos <= 33448354))
# 
# # key step. replace the original SNP column with the new column INCLUDING proxy
# exp_dat <- MR_LeadSNPs_proxy_5 %>% dplyr::select(-SNP) %>% dplyr::rename(SNP = SNP_out)   # finally, update exp_dat with the one with prox
# 
# } else { 
# #######################
# # in case of empty snappy results, no need to modify exp_dat. so do nothing
# #######################
#   print('no proxy was found by snappy')
# } # close if('proxy' %in% proxy_snp_cov_filtered_snp$match_check ){ 
# } # close if(file.exists(paste0('/scratch/richards/satoshi.yoshiji/18.obesity_cardiometab/00.fatpercentage_to_deCODE_proxy/snappy/', protname, '/snappyout.txt'))){

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
  saveRDS(mrpresso_res, file = 'output_COL6A3/mrpresso.RDS')
  
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
