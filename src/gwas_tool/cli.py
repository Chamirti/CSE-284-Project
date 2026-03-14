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
    
    parser.add_argument("--raw", required=True, help="Path to .raw genotype file")
    parser.add_argument("--bim", required=True, help="Path to .bim file for SNP positions")
    parser.add_argument("--causal", required=True, help="Path to causal.snplist for validation")
    parser.add_argument("--mode", choices=['naive', 'pca'], required=True, 
                        help="Choose 'naive' or 'pca' (includes LD pruning and ancestry correction)")
    
    args = parser.parse_args()

    display_name = "pca corrected with ld pruning" if args.mode == 'pca' else "naive"
    folder_name = display_name.replace(" ", "_")
    output_dir = os.path.join("results", folder_name)
    
    os.makedirs(output_dir, exist_ok=True)

    print(f"[*] Loading data from {args.raw}...")
    genotypes = pd.read_csv(args.raw, sep=r'\s+', low_memory=False)
    
    y = genotypes.iloc[:, 5].values
    X = genotypes.iloc[:, 6:].values
    
    col_means = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_means, inds[1])

    PCs = None
    if args.mode == 'pca':
        print(f"\n[*] Starting {display_name.upper()} pipeline...")
        print("[*] Step 1: Running LD Pruning (Window: 50, Step: 5, R² < 0.2)...")
        kept_indices = run_ld_pruning(X)
        
        print(f"[*] Step 2: Computing PCA via SVD on {len(kept_indices)} pruned SNPs...")
        PCs = run_pca(X[:, kept_indices], n_pcs=3)
        print("[*] PCA complete. Top 3 components extracted.")
    else:
        print(f"\n[*] Running {display_name.upper()} GWAS (No correction)...")

    betas, p_vals = run_gwas_math(X, y, PCs)

    snp_names = [name.split('_')[0] for name in genotypes.columns[6:]]
    results = pd.DataFrame({
        'SNP': snp_names, 
        'BETA': betas, 
        'P': p_vals
    }).dropna()
    
    chisq = stats.chi2.ppf(1 - results['P'], 1)
    lambda_gc = np.median(chisq) / 0.454

    causal_list = pd.read_csv(args.causal, header=None)[0].values
    false_positives = results[(~results['SNP'].isin(causal_list)) & (results['P'] < 5e-8)]

    print("\n" + "="*50)
    print(f"GWAS ANALYSIS REPORT: {display_name.upper()}")
    print("-" * 50)
    print(f"Genomic Inflation (λ_GC):   {lambda_gc:.3f}")
    print(f"False Positives Detected:   {len(false_positives)}")
    print(f"Output Directory:           {output_dir}")
    print("="*50)

    results.to_csv(f"{output_dir}/gwas_results.csv", index=False)
    generate_visuals(results, args.bim, causal_list, output_dir)
    
    print(f"\n[+] Success! Visualizations and CSV saved to {output_dir}\n")

if __name__ == "__main__":
    main()
