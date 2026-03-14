# CSE-284 Project  
# GWAS Tool: A Command-Line Genome-Wide Association Study (GWAS) Pipeline with LD Pruning and PCA Correction for Linear Phenotypes

![Python](https://img.shields.io/badge/python-3.10+-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Course](https://img.shields.io/badge/course-CSE284-orange)
![Status](https://img.shields.io/badge/status-course_project-lightgrey)

---

## Overview

This project implements a command-line Genome-Wide Association Study (GWAS) tool for analyzing linear (continuous) phenotypes. The tool performs association testing between genetic variants (SNPs) and a quantitative trait using Ordinary Least Squares (OLS) linear regression.

It supports two analysis modes:

- **Naive Mode** – Performs linear regression without population correction.  
- **PCA Mode** – Applies LD pruning followed by Principal Component Analysis (PCA) and includes the top principal components as covariates to correct for population stratification.

---

## Features

- Supports PLINK `.raw` genotype files  
- Processes phenotype data files  
- Automatic sample alignment using Individual IDs (IID)  
- Mean imputation for missing genotypes  
- Linkage Disequilibrium (LD) pruning  
- PCA-based population structure correction (top 3 PCs)  
- Linear regression association testing  
- Outputs SNP, Beta, p-value, and genomic inflation factor (λGC)  
- Generates Manhattan and Q–Q plots  
- Benchmark comparison with PLINK  

---

## Requirements

Dependencies are managed using `pyproject.toml`.

The project uses:

- numpy  
- pandas  
- scipy  
- matplotlib  
- seaborn  
- psutil  

No manual installation of dependencies is required if the project is installed using the provided build configuration.

---

## Quick Start (Example Usage)

### 1. Clone the Repository

```bash
git clone https://github.com/Chamirti/CSE-284-Project.git
```
```bash
cd CSE-284-Project
```

### 2. Install the Tool

Since this project uses `pyproject.toml`, install it with:

```bash
pip install -e .
```

### 3. Run the Tool

#### Naive GWAS:

```bash
python -m gwas_tool.cli \
  --raw example_data/genotypes.raw \
  --bim example_data/genotypes.bim \
  --causal example_data/causal.snplist \
  --mode naive
```

#### PCA Mode:

```bash
python -m gwas_tool.cli \
  --raw example_data/genotypes.raw \
  --bim example_data/genotypes.bim \
  --causal example_data/causal.snplist \
  --mode pca
```
Note: The provided example dataset is a small subset (100 individuals, 1000 SNPs) and is intended for testing purposes. It may not produce biologically meaningful associations.

---

## Outputs

Results are saved in:

- `Results/naive`
- `Results/pca_corrected_with_ld_pruning`

Each folder contains:

- Association statistics (SNP, Beta, p-value)
- Manhattan plot
- Q–Q plot


---

## Benchmarking

The tool was evaluated on a dataset containing:

- 50,000 SNPs  
- 310 individuals  

Results were benchmarked against **PLINK v1.90b** using:

```bash
--linear
--allow-no-sex
```

Performance was evaluated in terms of:

- Execution time  
- Memory usage  
- P-value correlation  
- Beta correlation  
- Genomic inflation factor (λGC)  

Benchmark implementation details are available in the **Benchmarking Results** folder.

---

## Implementation Details

- Python-based modular architecture  
- NumPy for vectorized matrix computations  
- Pandas for data integration  
- OLS regression implemented using linear algebra  
- LD pruning using a sliding window approach  
- PCA computed using Singular Value Decomposition (SVD)  

---

## Future Improvements

- Support for PLINK binary formats (.bed/.bim/.fam)  
- Additional covariates (age, sex, environmental factors)  
- Automated quality control (MAF, HWE filtering)  
- Chunked processing for large-scale datasets  
- Improved scalability for biobank-scale data  

---

## Author

Developed by Chamirti Senthilkumar as part of the CSE 284 course project.

---
