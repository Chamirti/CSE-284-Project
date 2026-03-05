#!/usr/bin/env bash
set -euo pipefail

# from repo root
python -m pip install -e .

gwas-tool gwas \
  --geno data/test/geno.csv \
  --pheno data/test/pheno.csv \
  --binary \
  --plots \
  --out results/test_results.csv

echo "Done. See results/test_results.csv and PNG plots."
