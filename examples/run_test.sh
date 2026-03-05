#!/bin/bash
set -e

echo "Running GWAS example (no PCA)..."

gwas-tool gwas \
  --geno data/test/geno.csv \
  --pheno data/test/pheno.csv \
  --binary \
  --plots \
  --out results/test_results.csv


echo "Running GWAS example (with PCA correction, 2 PCs)..."

gwas-tool gwas \
  --geno data/test/geno.csv \
  --pheno data/test/pheno.csv \
  --binary \
  --pcs 2 \
  --plots \
  --out results/test_results_pca.csv


echo "Running ancestry-stratified GWAS (no PCA)..."

gwas-tool gwas \
  --geno data/test/geno.csv \
  --pheno data/test/pheno.csv \
  --binary \
  --ancestry data/test/ancestry.csv \
  --out results/test_results_ancestry.csv


echo "Example complete. Check the results/ folder."
