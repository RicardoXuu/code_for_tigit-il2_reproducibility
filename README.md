# αTIGIT-IL2 single-cell analysis – code for reproducibility
> Reproduces the single-cell RNA-seq processing, trajectory analysis, differential
expression and visualisations used in  
> **“αTIGIT-IL2 achieves tumor regression by promoting tumor-infiltrating regulatory T-cell fragility”**  
> *(submmited)*

## Available analysis folders

| folder | main purpose | datasets analysed |
| ------ | ------------ | ----------------- |
| **analysis_for_original_data** | Reproduces every figure generated from the *new* single-cell RNA-seq produced for our study.  Pipeline steps mirror those described in the Methods section and use only the raw counts deposited in NGDC-GSA. | Raw FASTQ and matrix files: **CRA019549** (NGDC-GSA) |
| **analysis_for_public_data** | Benchmarks αTIGIT-IL2 signatures in previously published human & mouse immunotherapy datasets, and interrogates bulk cohorts for translational relevance.  Includes all QC, integration and figure scripts needed to regenerate Extended Data panels and Supplementary results. | - **GSE169246** (Zhang *et al.*, *Cancer Cell*, TNBC scRNA-seq & scATAC-seq)  
- **GSE136206** (Hollern *et al.* 2019, TNBC mouse scRNA-seq)  
- **GSE123814** (Yost *et al.* 2019, human basal-cell carcinoma scRNA-seq)  
- TCGA bulk RNA-seq for LUAD, SKCM, BRCA (downloaded from the **UCSC Xena** portal) |

> **Note** Raw human FASTQ files from Zhang *et al.* are privacy-restricted; all analyses here start from the processed matrices supplied by the authors.  
> All scripts can be executed without modification once the above data are placed in the expected sub-folders (see each script header for exact paths).
