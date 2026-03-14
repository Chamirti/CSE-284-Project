import argparse
import os
import pandas as pd
import numpy as np
from scipy import stats
from .gwas import run_gwas_math
from .pca_and_ld import run_ld_pruning, run_pca
from .plots import generate_visuals

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", required=True)
    parser.add_argument("--bim", required=True)
    parser.add_argument("--causal", required=True)
    parser.add_argument("--mode", choices=['naive', 'corrected'], required=True)
    args = parser.parse_args()

    # Create results folder
    os.makedirs("results", exist_ok=True)

    # Load Data
    genotypes = pd.read_csv(args.raw, sep=r'\s+', low_memory=False)
    y = genotypes.iloc[:, 5].values
    X = genotypes.iloc[:, 6:].values
    X[np.isnan(X)] = np.take(np.nanmean(X, axis=0), np.where(np.isnan(X))[1])

    PCs = None
    if args.mode == 'corrected':
        print("[*] Running LD Pruning and PCA...")
        kept = run_ld_pruning(X)
        PCs = run_pca(X[:, kept])

    print(f"[*] Running {args.mode} GWAS...")
    betas, p_vals = run_gwas_math(X, y, PCs)

    # Prepare Results
    snp_names = [n.split('_')[0] for n in genotypes.columns[6:]]
    results = pd.DataFrame({'SNP': snp_names, 'BETA': betas, 'P': p_vals}).dropna()
    
    # Validation Stats
    chisq = stats.chi2.ppf(1 - results['P'], 1)
    lambda_gc = np.median(chisq) / 0.454
    causal_list = pd.read_csv(args.causal, header=None)[0].values
    fps = len(results[(~results['SNP'].isin(causal_list)) & (results['P'] < 5e-8)])

    print("\n" + "="*30)
    print(f"Genomic Inflation (λ): {lambda_gc:.3f}")
    print(f"False Positives: {fps}")
    print("="*30)

    # Save outputs
    results.to_excel("results/gwas_results.xlsx", index=False)
    generate_visuals(results, args.bim, causal_list, "results")
    print("\nDone! Check the 'results' folder for your Excel and Plots.")

if __name__ == "__main__":
    main()
