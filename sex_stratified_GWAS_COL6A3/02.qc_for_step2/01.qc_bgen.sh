#!/bin/bash
#SBATCH --job-name=qc_bgen
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH -t 12:00:00
#SBATCH -o ./log/output.%j.out
#SBATCH -e ./log/output.%j.err
#SBATCH --array=5

cd $SLURM_SUBMIT_DIR
chr=${SLURM_ARRAY_TASK_ID}
mkdir -p per_chr_bgen_files

export PATH=/home/richards/satoshi.yoshiji/scratch/tools/plink2:$PATH

plink2 --bgen path-to-ukb-bgen/imputed.v3/bgen/${chr}.bgen 'ref-first' \
	   --sample path-to-ukb-bgen-sample/ukb27449_imp_chr1_v3_s487395_20210201.sample \
	   --autosome --mac 50 \
	   --extract maf0.0001_info0.8.snp.lst \
	   --export bgen-1.2 \
	   --out per_chr_bgen_files/${chr}.qc

