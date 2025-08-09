[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16785532.svg)](https://doi.org/10.5281/zenodo.16785532)
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

> **Note** Raw human FASTQ files from Zhang *et al.* are privacy-restricted; all analyses here start from the processed matrices supplied by the authors.  
