import argparse
import os
import pandas as pd
import numpy as np
from scipy import stats
from .gwas import run_gwas_math
from .pca_and_ld import run_ld_pruning, run_pca
from .plots import generate_visuals

def main():
    parser = argparse.ArgumentParser(description="gwas-tool: High-performance GWAS with PCA & LD Pruning")
    
    # Updated Arguments with Help text
    parser.add_argument("--raw", required=True, help="Path to the .raw genotype file (e.g., data/input.raw)")
    parser.add_argument("--bim", required=True, help="Path to the .bim file for SNP positions")
    parser.add_argument("--causal", required=True, help="Path to causal.snplist for validation")
    parser.add_argument("--mode", choices=['naive', 'corrected'], required=True, help="Choose 'naive' for simple GWAS or 'corrected' for PCA-based GWAS")
    
    args = parser.parse_args()

    # 1. Create results folder
    os.makedirs("results", exist_ok=True)

    # 2. Load Data from the user-provided paths
    genotypes = pd.read_csv(args.raw, sep=r'\s+', low_memory=False)
    y = genotypes.iloc[:, 5].values
    X = genotypes.iloc[:, 6:].values
    
    # Impute missing values
    X[np.isnan(X)] = np.take(np.nanmean(X, axis=0), np.where(np.isnan(X))[1])

    PCs = None
    if args.mode == 'corrected':
        print("[*] Running LD Pruning and PCA calculation...")
        kept = run_ld_pruning(X)
        PCs = run_pca(X[:, kept])

    print(f"[*] Running {args.mode} GWAS math...")
    betas, p_vals = run_gwas_math(X, y, PCs)

    # 3. Process Results
    snp_names = [n.split('_')[0] for n in genotypes.columns[6:]]
    results = pd.DataFrame({'SNP': snp_names, 'BETA': betas, 'P': p_vals}).dropna()
    
    # 4. Validation Stats (Printed to CMD)
    chisq = stats.chi2.ppf(1 - results['P'], 1)
    lambda_gc = np.median(chisq) / 0.454
    
    causal_list = pd.read_csv(args.causal, header=None)[0].values
    fps = len(results[(~results['SNP'].isin(causal_list)) & (results['P'] < 5e-8)])

    print("\n" + "="*40)
    print(f"VALDIATION METRICS")
    print("-" * 40)
    print(f"Genomic Inflation (λ_GC): {lambda_gc:.3f}")
    print(f"False Positives Detected: {fps}")
    print("="*40)

    # 5. Save Outputs to the 'results' folder
    results.to_csv("results/gwas_results.csv", index=False)
    generate_visuals(results, args.bim, causal_list, "results")
    
    print("\n[+] Analysis Complete!")
    print("Results saved: results/gwas_results.csv")
    print("Visuals saved: results/gwas_visualizations.png")

if __name__ == "__main__":
    main()
