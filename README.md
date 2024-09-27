---
## COL6A3-Derived Endotrophin Mediates the Effect of Obesity on Coronary Artery Disease: An Integrative Proteogenomics Analysis  
*medRxiv 2023*  
[DOI: https://doi.org/10.1101/2023.04.19.23288706](https://doi.org/10.1101/2023.04.19.23288706)

---

### Step 1: MR - Identifying the Effect of BMI on Plasma Protein Levels  
1. Generate instrumental variables by running `00.remove_MHC_and_clump.R` in the `1.1.BMI_to_proteins/` directory.  
2. Run `01.MR.R`.

### Step 2: MR - Identifying the Effect of BMI-Driven Proteins on Cardiometabolic Outcomes  
1. Run `01.MR.R` in the `2.1.protein_to_disease` directory.

### Colocalization Analysis  
1. Run `01.coloc.R` in the `2.2.colocalization` directory.  
2. `grid.csv` provides parameters for prior probabilities used in the coloc analysis.

### Sex-Stratified GWAS of C-terminal COL6A3 Levels in the UK Biobank
1. Run `01.filter.sh` and `02.prune.sh` in `02.qc_for_step1` to perform QC on variants for regenie step 1 (using merged unimputed variants).
2. Run `01.qc_bgen.sh` in `02.qc_for_step2.sh` to QC variants per chromosome for regenie step 2 (using imputed variants).
3. Run `01.regenie_step1.sh` in each directory within `03.step1` to execute regenie step 1.
4. Run `01.regenie_step2.sh` in each directory within `04.step2` to execute regenie step 2.
