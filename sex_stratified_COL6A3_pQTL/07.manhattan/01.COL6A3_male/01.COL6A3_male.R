library(tidyverse)
library(data.table)
#library(qqman)
library(fastman)
library(fastqq)

##########
# other options
# https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
##########

setwd("/home/richards/satoshi.yoshiji/scratch/30.pQTL_regenie/07.manhattan/01.COL6A3_male/")

gwas_all <- data.frame()
for(i in 1:22){
  gwas_chr <- fread(paste0('/scratch/richards/satoshi.yoshiji/30.pQTL_regenie/06.step2/01.COL6A3_male/output/COL6A3.chr', i, '_Protein_level.regenie'))
  gwas_all <- rbind(gwas_all, gwas_chr)
}

gwas_all_p <- gwas_all %>% mutate(pval = 10^(-LOG10P))

# check
# min(gwas_all_p$pval)
# max(gwas_all_p$pval)

# remove variants with pval >= 0.05 to ease plotting
gwas_cut <- gwas_all_p %>% filter(pval < 0.05)
write_tsv(gwas_cut, file = 'gwas_cut.tsv.gz')
write_tsv(gwas_all_p, file = 'gwas_all.tsv.gz') # also export the full gwas¬

# GC lambda
(median_pval <- median(gwas_all_p$pval))
chisq <- qchisq(median_pval, 1, lower.tail = F)
(gclambda <- chisq/0.456)

write.table(gclambda, file = 'gclambda.txt', quote = F, row.names = F, col.names = F)
  
# manhattan
pdf('manhattan.pdf', width = 8, height = 4)
fastman(gwas_cut, chr = "CHROM", bp = "GENPOS", p = "LOG10P", snp = 'ID', 
        logp = F, maxP = NULL, speedup = T, 
        suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8), baseline = 0,
        ylim = c(0, 35))
dev.off()

# pdf('manhattan.pdf', width = 12, height = 6)
# manhattan(gwas_cut, chr = 'CHROM', bp = 'GENPOS', p = 'pval', snp = 'ID')
# dev.off()

pdf('qq.pdf', width = 6, height = 6)
fastqq::qq(gwas_all_p$pval, main = "Q-Q plot of GWAS p-values")
dev.off()
