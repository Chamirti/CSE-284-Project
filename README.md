# CSE-284-Project  
# GWAS Tool (Course Project)

![Python](https://img.shields.io/badge/python-3.10+-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Course](https://img.shields.io/badge/course-CSE284-orange)
![Status](https://img.shields.io/badge/status-course_project-lightgrey)

Lightweight command-line GWAS tool supporting PCA correction, ancestry-stratified analysis, and Manhattan/QQ plots.
This repository contains a lightweight **command-line GWAS (Genome-Wide Association Study) implementation** developed for the CSE 284 course project.

The tool performs **SNP-wise association testing** using genotype and phenotype data provided as CSV files. It also includes optional **population structure correction using PCA**, **ancestry-stratified GWAS**, and visualization of results using **Manhattan and QQ plots**.

The project is designed as an **educational implementation** to demonstrate the core components of a GWAS pipeline on a small toy dataset.

---

# Pipeline Overview

The tool implements a simplified GWAS pipeline consisting of the following steps:

1. **Load input data**
   - Genotype matrix (`samples × SNPs`) from CSV  
   - Phenotype vector (`sample_id`, `phenotype`)

2. **Align samples**  
   Genotype and phenotype files are merged using the `sample_id` column.

3. **Optional PCA correction**  
   Principal components are computed from the genotype matrix and included as covariates in the regression model to correct for population structure.

4. **SNP-wise association testing**  
   Each SNP is tested independently using regression:

   - **Logistic regression** for binary phenotypes  
   - **Linear regression** for continuous phenotypes

5. **Result generation**  
   The tool outputs a table containing:
   - SNP identifier  
   - effect size (beta)  
   - p-value

6. **Visualization**
   - Manhattan plot
   - QQ plot

7. **Optional ancestry-stratified GWAS**  
   If an ancestry file is provided, the dataset is split by ancestry group and GWAS is run separately for each group.

Pipeline summary:

    Genotype + Phenotype
            │
            ▼
       Align samples
            │
            ▼
       Optional PCA
            │
            ▼
       SNP-wise GWAS
            │
            ├── Manhattan plot
            └── QQ plot
            │
            ▼
    Ancestry-stratified GWAS (optional)

---

# Features

### GWAS analysis
- SNP-wise association testing
- Supports:
  - Linear regression (continuous phenotypes)
  - Logistic regression (binary phenotypes)

Outputs:
- SNP identifier
- Effect size (beta)
- p-value

### Population structure correction
- Optional **PCA covariates**
- Principal components computed directly from genotype matrix
- Helps control for population stratification

### Ancestry-stratified GWAS
- Optional ancestry file (`sample_id`, `ancestry`)
- Runs separate GWAS analyses for each ancestry group
- Outputs separate result files for each ancestry

### Visualization
Automatically generates:

- Manhattan plot  
- QQ plot

(using SNP index for the toy dataset)

### Example dataset
Repository includes a small **toy dataset (10 samples, 5 SNPs)** to demonstrate the pipeline.

---

# Installation

Install the package in editable mode:

    python -m pip install -e .

This installs the command-line tool:

    gwas-tool --help

---
# Quick Start (Example Usage):

### 1. Basic GWAS

    gwas-tool gwas \
      --geno data/test/geno.csv \
      --pheno data/test/pheno.csv \
      --binary \
      --plots \
      --out results/test_results.csv

Runs SNP-wise association testing and generates Manhattan and QQ plots.

---

### 2. PCA-Corrected GWAS

    gwas-tool gwas \
      --geno data/test/geno.csv \
      --pheno data/test/pheno.csv \
      --binary \
      --pcs 2 \
      --plots \
      --out results/test_results_pca.csv

Adds the top **2 principal components** as covariates to correct for population structure.

---

### 3. Ancestry-Stratified GWAS

    gwas-tool gwas \
      --geno data/test/geno.csv \
      --pheno data/test/pheno.csv \
      --binary \
      --ancestry data/test/ancestry.csv \
      --out results/test_results_ancestry.csv

Runs GWAS separately for each ancestry group and outputs ancestry-specific results files.

This produces separate results files for each ancestry group.

Example outputs:

    results/
      test_results_ancestry_AFR.csv
      test_results_ancestry_EUR.csv
      test_results_ancestry_EAS.csv

---

# Running the Example Script

The repository includes a demonstration script:

    examples/run_test.sh

Run:

    bash examples/run_test.sh

This script runs:

- Standard GWAS
- PCA-corrected GWAS
- Ancestry-stratified GWAS

and generates results in the `results/` directory.

---

# Example Output

GWAS results file:

    snp,beta,p_value
    SNP_A,0.29,0.57
    SNP_B,-0.29,0.57
    ...

Plots generated:

    results/test_results.manhattan.png
    results/test_results.qq.png

---

# Project Structure

    CSE-284-Project
    │
    ├── src/gwas_tool
    │   ├── cli.py
    │   ├── gwas.py
    │   ├── pca.py
    │   └── plots.py
    │
    ├── data/test
    │   ├── geno.csv
    │   ├── pheno.csv
    │   └── ancestry.csv
    │
    ├── examples
    │   └── run_test.sh
    │
    ├── pyproject.toml
    └── README.md

---

# Limitations

- Example dataset is intentionally **very small** (10 samples, 5 SNPs)
- P-values may appear unstable for ancestry-stratified analyses due to small sample sizes
- Intended for **educational demonstration**, not large-scale GWAS analysis

---

# Remaining Work

Planned extension:

- Linkage Disequilibrium (LD) analysis / SNP correlation exploration

---
