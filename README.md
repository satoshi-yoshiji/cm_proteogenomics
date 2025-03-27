## Integrative proteogenomic analysis identifies COL6A3-derived endotrophin as a mediator of the effect of obesity on coronary artery disease  
Satoshi Yoshiji, Tianyuan Lu, Guillaume Butler-Laporte, Julia Carrasco-Zanini-Sanchez, Chen-Yang Su, Yiheng Chen, Kevin Liang, Julian Daniel Sunday Willett, Shidong Wang, Darin Adra, Yann Ilboudo, Takayoshi Sasako, Satoshi Koyama, Tetsushi Nakao, Vincenzo Forgetta, Yossi Farjoun, Hugo Zeberg, Sirui Zhou, Michael Marks-Hultström, Mitchell J. Machiela, Rama Kaalia, Hesam Dashti, Melina Claussnitzer, Jason Flannick, Nicholas J. Wareham, Vincent Mooser, Nicholas J. Timpson, Claudia Langenberg & J. Brent Richards.  
*Nature Genetics* **57**, 345–357 (2025).  
[https://doi.org/10.1038/s41588-024-02052-7](https://doi.org/10.1038/s41588-024-02052-7)


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
