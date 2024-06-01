#!/bin/bash
#SBATCH --job-name=step2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96GB
#SBATCH -t 36:00:00
#SBATCH -o ./log/output.%j.out
#SBATCH -e ./log/output.%j.err
#SBATCH --array 2
cd $SLURM_SUBMIT_DIR

# conda
source ~/anaconda3/etc/profile.d/conda.sh
conda activate regenie

mkdir -p output/

regenie \
  --step 2 \
  --bgen /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/04.merge_imputed_files/per_chr_bgen_files/${SLURM_ARRAY_TASK_ID}.qc.bgen \
  --sample /project/richards/restricted/ukb-27449/data/ukb27449_imp_chr1_v3_s487395_20210201.sample \
  --keep /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_female_RNT/02.keep_female.tsv \
  --chr ${SLURM_ARRAY_TASK_ID} \
  --covarFile /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_female_RNT/02.covariate_female.tsv \
  --phenoFile /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/03.pheno/01.COL6A3_female_RNT/02.pheno_female.tsv \
  --range 2:238232752-238232752 \
  --bsize 200 \
  --phenoCol Protein_level \
  --threads 20 \
  --pred /scratch/richards/satoshi.yoshiji/30.pQTL_regenie/05.step1/01.COL6A3_female_RNT/fit_bin_out_pred.list --ref-first \
  --out output/COL6A3.chr${SLURM_ARRAY_TASK_ID}
  
# regenie does not accept both --keep and --remove. Thus, I used R to exclude some individuals from the keep list
#--remove /scratch/richards/satoshi.yoshiji/08.UKB/excluded_individuals/w27449_2023-04-25.tsv \
