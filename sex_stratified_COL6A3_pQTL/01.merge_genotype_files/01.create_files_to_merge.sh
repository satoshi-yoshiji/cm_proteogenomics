#!/bin/bash

mkdir -p per_chr_plink_files

# copy bed files (symbolic link)_
for i in {1..22};do ln -s /scratch/richards/restricted/ukb-general/scratch/genetic-data/genome/bed/${i}.bed  per_chr_plink_files/${i}.bed ; done

# copy bim files (symbolic link)_
for i in {1..22};do ln -s /scratch/richards/restricted/ukb-general/scratch/genetic-data/genome/bed/${i}.bim  per_chr_plink_files/${i}.bim ; done

# copy fam files (the same for all 22 chromosomes)
for i in {1..22};do cp /scratch/richards/restricted/ukb-27449/scratch/old-storage/full_release/v3/genotyped/w27449_20200204/ukb27449_cal_chr1_v2_20200204.fixCol6.fam  per_chr_plink_files/${i}.fam ; done

# create a list of files to be merged
ls per_chr_plink_files/[0-9]*.bed | sort -V | sed -e 's/\.bed//g' > files_to_merge.txt
