# colocalization for step 2

# Load necessary libraries
library(vroom)
library(coloc)
library(tidyverse)
library(magrittr)

# working directory
wd <- 'path-to-wd'
setwd(wd)

# Create output directory for results
system('mkdir -p 04.pwcoco/pwcoco_gwas/')

# Parse command line arguments: outcome_name, outcome_path, protname, samplesize, casesize
args <- commandArgs(trailingOnly = TRUE)
outcome_name <- args[1]
outcome_path <- args[2]
protname <- args[3]
samplesize <- as.numeric(args[4])
casesize <- as.numeric(args[5])

# Print sample size and case size for verification
print(samplesize)
print(casesize)

# Extract sequence ID from protein name
seqid <- str_split(protname, '[.]')[[1]][2]

# Load cis-pQTL data
cis_pqtl_path <- system(paste0('ls path-to-cis-pQTL/*', seqid, '*.tsv'), intern = TRUE)

if (file.exists(cis_pqtl_path)) {
  tryCatch({
    # Read cis-pQTL file and get the top variant based on p-value
    cispqtl <- vroom(cis_pqtl_path) %>%
      arrange(desc(log10p_unadj)) %>%
      slice(1)
    
    # Extract relevant information (rsid, chr, position)
    cispqtl_lowestpval_rsid <- cispqtl %>% select(variant)
    cispqtl_lowestpval_chr <- cispqtl %>% select(chr_var) %>% unlist()
    cispqtl_lowestpval_pos <- cispqtl %>% select(pos_var) %>% unlist()
    cispqtl_lowestpval_chr_num <- cispqtl %>%
      select(chr_var) %>%
      separate(chr_var, into = c('temp', 'chr'), sep = 'r') %>%
      select(chr) %>%
      as.numeric()
    
    # Load corresponding pGWAS data
    pgwas_path <- system(paste0('ls path-to-pQTL/*', seqid, '*.txt.gz'), intern = TRUE)
    pgwas <- vroom(pgwas_path) %>%
      filter(Chrom == cispqtl_lowestpval_chr) %>%
      filter(Pos >= cispqtl_lowestpval_pos - 250000, Pos <= cispqtl_lowestpval_pos + 250000) %>%
      mutate(varbeta = SE^2, snp = unlist(lapply(strsplit(rsids, ","), function(x) x[1]))) %>%
      filter(!is.na(snp)) %>%
      arrange(Pval) %>%
      filter(!duplicated(snp))
    
    # Rename columns in pGWAS
    pgwas <- pgwas %>%
      rename(
        position = Pos,
        beta = Beta,
        varbeta = varbeta,
        ALT = effectAllele,
        REF = otherAllele,
        MAF = ImpMAF,
        pvalues = Pval
      ) %>%
      separate(Chrom, into = c('former', 'chr'), sep = 'r') %>%
      select(-former) %>%
      unite(col = 'id', chr, position, REF, ALT, sep = ":", remove = FALSE) %>%
      drop_na()
    
    # Handle edge cases for MAF and p-value
    pgwas <- pgwas %>%
      mutate(
        MAF = ifelse(MAF == 0, 1e-4, ifelse(MAF == 1, 0.9999, MAF)),
        pvalues = ifelse(pvalues == 0, 1e-300, ifelse(pvalues == 1, 0.9999, pvalues))
      )
    
    # Load outcome GWAS data and filter for relevant SNPs
    outcome <- vroom(outcome_path) %>%
      rename(beta = ES, SE = SE, pvalues = pval, MAF = AF, snp = ID) %>%
      filter(snp %in% unlist(pgwas$snp)) %>%
      mutate(varbeta = SE^2) %>%
      arrange(pvalues) %>%
      filter(!duplicated(snp)) %>%
      mutate(
        MAF = ifelse(MAF == 0, 1e-4, ifelse(MAF == 1, 0.9999, MAF)),
        pvalues = ifelse(pvalues == 0, 1e-300, ifelse(pvalues == 1, 0.9999, pvalues))
      )
    
    # Identify common SNPs between pGWAS and outcome data
    common <- inner_join(pgwas, outcome, by = 'snp', suffix = c('.pgwas', '.outcome'))
    
    if (nrow(common) == 0) {
      # If no common SNPs, write an error message and stop the process
      resname <- paste0(protname, '.', outcome_name, '.coloc.tsv')
      write.table('no common snp', file = resname, quote = FALSE, row.names = FALSE, sep = '\t')
      print('no common snp')
    } else {
      # Save common SNP results
      pgwas_common_filename <- paste0('04.pwcoco/pwcoco_gwas/', protname, '.', outcome_name, '.common.500kb.tsv')
      outcome_common_filename <- paste0('04.pwcoco/pwcoco_gwas/', outcome_name, '.', protname, '.common.500kb.tsv')
      write.table(common %>% select(snp, ALT.pgwas, REF.pgwas, beta.pgwas, SE.pgwas, pvalues.pgwas, MAF.pgwas), file = pgwas_common_filename, sep = '\t', quote = FALSE, row.names = FALSE)
      write.table(common %>% select(snp, ALT.outcome, REF.outcome, beta.outcome, SE.outcome, pvalues.outcome, MAF.outcome), file = outcome_common_filename, sep = '\t', quote = FALSE, row.names = FALSE)
      
      # Prepare data for coloc
      sumstats1_list <- list(
        type = 'quant',
        N = 35559,
        MAF = common$MAF.pgwas,
        beta = common$beta.pgwas,
        varbeta = common$varbeta.pgwas,
        pvalues = common$pvalues.pgwas
      )
      
      sumstats2_list <- list(
        type = 'cc',
        N = samplesize,
        s = casesize / samplesize,
        MAF = common$MAF.outcome,
        beta = common$beta.outcome,
        varbeta = common$varbeta.outcome,
        pvalues = common$pvalues.outcome
      )
      
      # Perform colocalization analysis using coloc.abf
      coloc_res <- coloc.abf(dataset1 = sumstats1_list, dataset2 = sumstats2_list)
      
      # Save colocalization results
      system('mkdir -p 04.coloc_res/')
      resname <- paste0('04.coloc_res/', protname, '-', outcome_name, '-coloc.tsv')
      write_tsv(as.data.frame(coloc_res$summary) %>% rownames_to_column(), file = resname)
    }
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
} else {
  print('No cis-pQTL file found.')
}
