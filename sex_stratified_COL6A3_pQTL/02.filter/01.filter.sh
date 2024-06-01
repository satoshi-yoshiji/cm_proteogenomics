#!/bin/bash

plink2 --bfile /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/01.merge_genotype_files/merged_plink_files/ukb_genotype_merged \
		   --autosome --maf 0.01 --geno 0.01 --mac 100 --hwe 1e-15 --make-bed \
		   --out 01.filtered/filtered
