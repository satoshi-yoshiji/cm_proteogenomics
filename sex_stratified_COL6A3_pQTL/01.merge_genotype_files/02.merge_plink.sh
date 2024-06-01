#!/bin/bash

mkdir -p merged_plink_files

plink --merge-list files_to_merge.txt --make-bed --autosome --out merged_plink_files/ukb_genotype_merged
