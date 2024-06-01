#!/bin/bash
#SBATCH --job-name=step1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96GB
#SBATCH -t 24:00:00
#SBATCH -o ./log/output.%j.out
#SBATCH -e ./log/output.%j.err

cd $SLURM_SUBMIT_DIR

# conda
source ~/anaconda3/etc/profile.d/conda.sh
conda activate regenie

regenie \
  --step 1 \
  --bed /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/01.merge_genotype_files/merged_plink_files/ukb_genotype_merged \
  --extract /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/02.filter/02.pruned/pruned.prune.in \
  --keep /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_female_RNT/02.keep_female.tsv \
  --covarFile /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_female_RNT/02.covariate_female.tsv \
  --phenoFile /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_female_RNT/02.pheno_female.tsv \
  --bsize 1000 \
  --phenoCol Protein_level \
  --threads 16 \
  --out fit_bin_out
  
#--exclude example/snplist_rm.txt \
  #--remove example/fid_iid_to_remove.txt \

##############
# one liner to be run on the rstudio node
##############
# regenie --step 1 --bed /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/01.merge_genotype_files/merged_plink_files/ukb_genotype_merged --extract /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/02.filter/02.pruned/pruned.prune.in --phenoFile /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3/02.pheno.tsv --covarFile /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3/02.covariate.tsv --bsize 100 --out ukb_training --threads 16 --gz
