#!/bin/bash

plink2 --bfile ukb_genotype_merged \
		   --autosome --maf 0.01 --geno 0.01 --mac 100 --hwe 1e-15 --make-bed \
		   --out 01.filtered/filtered
