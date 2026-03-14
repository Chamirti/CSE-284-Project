# CSE-284-Project  
# GWAS Tool: A Command-Line Genome-Wide Association Study (GWAS) Pipeline with LD Pruning and PCA Correction for Linear Phenotypes

![Python](https://img.shields.io/badge/python-3.10+-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Course](https://img.shields.io/badge/course-CSE284-orange)
![Status](https://img.shields.io/badge/status-course_project-lightgrey)

## Overview

This project implements a Genome-Wide Association Study (GWAS) tool for analyzing linear phenotypes. The tool performs association testing between genetic variants (SNPs) and a continuous trait using Ordinary Least Squares (OLS) linear regression.

It includes two analysis modes:

- **Naive Mode** – Performs linear regression without population correction.
- **PCA Mode** – Applies LD pruning followed by Principal Component Analysis (PCA) and includes the top principal components as covariates to correct for population stratification.

---

## Features

- Supports PLINK `.raw` genotype files
- Handles phenotype data from `.phen` files
- Automatic sample alignment using Individual IDs (IID)
- Mean imputation for missing genotypes
- LD pruning before PCA
- Population structure correction using PCA (top 3 PCs)
- Linear regression association testing
- Outputs Beta, p-value, genomic inflation factor
- Generates Manhattan and Q–Q plots
- Benchmark comparison with PLINK

---

## Requirements

- Python 3.x
- NumPy (v1.x)
- Pandas (v2.x)
- Matplotlib
- SciPy
- scikit-learn

---

## Quick Start (Example Usage)

### 1. Clone the Repository

```bash
git clone <your-repository-link>
cd <your-project-folder>
```

### 2. Install Dependencies

```bash
pip install -r requirements.txt
```

(If no `requirements.txt` is provided, install dependencies manually using the command above.)

### 3. Prepare Input Files

You need:

- Genotype file in **PLINK `.raw` format**
- Phenotype file in **`.phen` format**

Ensure both files share matching **Individual IDs (IID)**.

### 4. Run the Tool

Example usage:

```bash
python gwas.py --genotype data.raw --phenotype data.phen --mode pca
```

Available modes:

- `--mode naive`
- `--mode pca`

### 5. View Results

The tool outputs:

- Association statistics (Beta, SE, t-statistic, p-value)
- Manhattan plot
- Q–Q plot
- Benchmark comparison metrics (if enabled)

---

## Benchmarking

The tool was evaluated on a dataset containing:

- 50,000 SNPs
- 310 individuals

Results were compared against **PLINK v1.90b** using:

```bash
--linear
--allow-no-sex
```

Performance was measured using:

- Execution time
- Memory usage
- P-value correlation
- Beta correlation
- Genomic inflation factor (λGC)

---

## Implementation Details

- Written in Python using a modular design
- Uses NumPy for vectorized matrix operations
- Uses Pandas for data integration
- OLS regression implemented using linear algebra
- LD pruning performed using sliding window approach
- PCA computed using Singular Value Decomposition (SVD)

---

## Future Improvements

- Support for PLINK binary formats (.bed/.bim/.fam)
- Additional covariates (age, sex, environmental factors)
- Automated quality control (MAF, HWE filtering)
- Chunked processing for large-scale datasets
- Enhanced scalability for biobank-sized data

---

## Author

Developed as part of a GWAS implementation project focusing on linear phenotype analysis with population structure correction.

---
