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

   # 1. Load Data
    print(f"[*] Loading genotypes from {args.raw}...")
    genotypes = pd.read_csv(args.raw, sep=r'\s+', low_memory=False)

    # NEW: Load the external binary phenotype
    # Assuming the file is named binary_gwas_data.phen and is in the same folder
    pheno_path = "C:/Users/chami/Downloads/binary_gwas_data.phen" 
    print(f"[*] Loading binary phenotypes from {pheno_path}...")
    
    # PLINK .phen files usually have FID, IID, and then the Trait
    pheno_df = pd.read_csv(pheno_path, sep=r'\s+', header=None, names=['FID', 'IID', 'y_binary'])

    # Merge on IID to ensure the DNA matches the correct Person
    merged_data = pd.merge(genotypes, pheno_df, on='IID')

    # Now extract y and X from the merged dataframe
    y = merged_data['y_binary'].values
    # Genotypes start after the metadata columns (FID, IID, PAT, MAT, SEX, PHENO)
    # Since we merged, check your column indices, but usually:
    X = merged_data.iloc[:, 7:-1].values # Adjust indices based on your merged df structure
    # 2. Simple Imputation
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
        print(f"\n[*] Running {display_name.upper()} GWAS...")

    # 3. Core Math Execution (Auto-switches between Linear and Logistic)
    betas, p_vals = run_gwas_math(X, y, PCs)

    # 4. Prepare Results DataFrame
    snp_names = [name.split('_')[0] for name in genotypes.columns[6:]]
    results = pd.DataFrame({
        'SNP': snp_names, 
        'BETA': betas, 
        'P': p_vals
    }).dropna()
    
    # 5. Calculate Validation Metrics (Lambda GC and False Positives)
    # Using Chi-squared for p-value inflation check
    chisq = stats.chi2.ppf(1 - results['P'], 1)
    lambda_gc = np.median(chisq) / 0.454
    
    # Load causal list to check accuracy
    causal_list = pd.read_csv(args.causal, header=None)[0].values
    false_positives = results[(~results['SNP'].isin(causal_list)) & (results['P'] < 5e-8)]

    # 6. Final Terminal Report
    print("\n" + "="*50)
    print(f"GWAS ANALYSIS REPORT: {display_name.upper()}")
    print("-" * 50)
    print(f"Trait Type:               {trait_type.upper()}")
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
