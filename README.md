# CSE-284-Project
# GWAS Tool (Course Project)

This repository contains a small command-line GWAS implementation for a course project.  
Current functionality supports SNP-wise association testing on CSV genotype/phenotype inputs and generates basic Manhattan and QQ plots.

## Features (current)
- SNP-wise GWAS using:
  - Linear regression (continuous phenotype)
  - Logistic regression (binary phenotype)
- Outputs per-SNP effect size (beta) and p-value
- Manhattan + QQ plots (toy: SNP index)

## Install
```bash
python -m pip install -e .
