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
  --bed path-to-wd/01.merge_genotype_files/merged_plink_files/ukb_genotype_merged \
  --extract path-to-wd/02.filter/02.pruned/pruned.prune.in \
  --keep path-to-wd/03.pheno/01.COL6A3_male_RNT/02.keep_male.tsv \
  --covarFile path-to-wd/03.pheno/01.COL6A3_male_RNT/02.covariate_male.tsv \
  --phenoFile path-to-wd/03.pheno/01.COL6A3_male_RNT/02.pheno_male.tsv \
  --bsize 1000 \
  --phenoCol Protein_level \
  --threads 16 \
  --out fit_bin_out