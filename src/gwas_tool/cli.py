import argparse
import os
import pandas as pd
import numpy as np
from scipy import stats

# Internal imports from your package modules
from .gwas import run_gwas_math
from .pca_and_ld import run_ld_pruning, run_pca
from .plots import generate_visuals

def main():
    parser = argparse.ArgumentParser(description="gwas-tool: High-performance GWAS with PCA & LD Pruning")
    
    # Arguments
    parser.add_argument("--raw", required=True, help="Path to .raw genotype file")
    parser.add_argument("--bim", required=True, help="Path to .bim file for SNP positions")
    parser.add_argument("--causal", required=True, help="Path to causal.snplist for validation")
    parser.add_argument("--mode", choices=['naive', 'pca'], required=True, 
                        help="Choose 'naive' or 'pca' (includes LD pruning and ancestry correction)")
    
    args = parser.parse_args()

    # 1. Define folder names and display labels
    display_name = "pca corrected with ld pruning" if args.mode == 'pca' else "naive"
    folder_name = display_name.replace(" ", "_")
    output_dir = os.path.join("results", folder_name)
    
    # Create the mode-specific results folder
    os.makedirs(output_dir, exist_ok=True)

    # 2. Load Data from the user-provided paths
    print(f"[*] Loading data from {args.raw}...")
    genotypes = pd.read_csv(args.raw, sep=r'\s+', low_memory=False)
    
    # Extract phenotype (y) and genotype matrix (X)
    y = genotypes.iloc[:, 5].values
    X = genotypes.iloc[:, 6:].values
    
    # Simple Imputation: Replace NaNs with the mean of the column
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

    # 3. Core Math Execution
    betas, p_vals = run_gwas_math(X, y, PCs)

    # 4. Prepare Results DataFrame
    snp_names = [name.split('_')[0] for name in genotypes.columns[6:]]
    results = pd.DataFrame({
        'SNP': snp_names, 
        'BETA': betas, 
        'P': p_vals
    }).dropna()
    
    # 5. Calculate Validation Metrics (Lambda GC and False Positives)
    chisq = stats.chi2.ppf(1 - results['P'], 1)
    lambda_gc = np.median(chisq) / 0.454
    
    # Load causal list to check accuracy
    causal_list = pd.read_csv(args.causal, header=None)[0].values
    false_positives = results[(~results['SNP'].isin(causal_list)) & (results['P'] < 5e-8)]

    # 6. Final Terminal Report
    print("\n" + "="*50)
    print(f"GWAS ANALYSIS REPORT: {display_name.upper()}")
    print("-" * 50)
    print(f"Genomic Inflation (λ_GC):   {lambda_gc:.3f}")
    print(f"False Positives Detected:   {len(false_positives)}")
    print(f"Output Directory:           {output_dir}")
    print("="*50)

    # 7. Save Outputs
    results.to_csv(f"{output_dir}/gwas_results.csv", index=False)
    generate_visuals(results, args.bim, causal_list, output_dir)
    
    print(f"\n[+] Success! Visualizations and CSV saved to {output_dir}\n")

if __name__ == "__main__":
    main()
