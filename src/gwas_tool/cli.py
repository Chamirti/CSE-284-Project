import argparse
from pathlib import Path
import pandas as pd

from .gwas import run_gwas_linear, run_gwas_logistic
from .pca import compute_pcs
from .plots import manhattan_plot, qq_plot


def main():
    p = argparse.ArgumentParser(prog="gwas-tool")
    sub = p.add_subparsers(dest="cmd", required=True)

    # GWAS command
    g = sub.add_parser("gwas", help="Run SNP-wise GWAS on geno/pheno CSV files.")
    g.add_argument("--geno", required=True, help="Genotype CSV: rows=samples, cols=SNPs (0/1/2), with first column sample_id.")
    g.add_argument("--pheno", required=True, help="Phenotype CSV: columns sample_id, phenotype (0/1 for binary or float for continuous).")
    g.add_argument("--out", required=True, help="Output results CSV path.")
    g.add_argument("--binary", action="store_true", help="Use logistic regression (phenotype must be 0/1). Default is linear regression.")
    g.add_argument("--pcs", type=int, default=0, help="Number of PCA covariates to include (0 = none).")
    g.add_argument("--plots", action="store_true", help="Create Manhattan and QQ plots next to output.")
    g.add_argument("--pval-col", default="p_value", help="Name of p-value column in output (default: p_value).")
    g.add_argument("--ancestry", help="Optional CSV file with columns: sample_id, ancestry")
    
    args = p.parse_args()

    geno_path = Path(args.geno)
    pheno_path = Path(args.pheno)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    geno = pd.read_csv(geno_path)
    pheno = pd.read_csv(pheno_path)
    ancestry = None
    if args.ancestry:
        ancestry = pd.read_csv(args.ancestry)

    # Align samples
    merged = pheno.merge(geno, on="sample_id", how="inner")
    if merged.shape[0] < 3:
        raise SystemExit("Not enough overlapping samples between geno and pheno.")

    y = merged["phenotype"].to_numpy()
    X = merged.drop(columns=["sample_id", "phenotype"])
    snp_names = list(X.columns)

    cov = None
    if args.pcs and args.pcs > 0:
        pcs_df = compute_pcs(X, n_components=args.pcs)
        cov = pcs_df.to_numpy()

    if args.binary:
        res = run_gwas_logistic(X.to_numpy(), y, snp_names, covariates=cov)
    else:
        res = run_gwas_linear(X.to_numpy(), y, snp_names, covariates=cov)

    res.to_csv(out_path, index=False)

    if args.plots:
        manhattan_plot(res, out_path.with_suffix(".manhattan.png"), pval_col=args.pval_col)
        qq_plot(res, out_path.with_suffix(".qq.png"), pval_col=args.pval_col)

    print(f"Wrote results: {out_path}")
    if args.plots:
        print(f"Wrote plots: {out_path.with_suffix('.manhattan.png')} and {out_path.with_suffix('.qq.png')}")

    ##ancestry stratified gwas
    if ancestry is not None:

        merged = pheno.merge(geno, on="sample_id").merge(ancestry, on="sample_id")
    
        for group in merged["ancestry"].unique():
    
            subset = merged[merged["ancestry"] == group]
    
            y = subset["phenotype"].to_numpy()

            if len(set(y)) < 2:
                print(f"Skipping ancestry group '{group}' because phenotype contains only one class.")
                continue
            
            X = subset.drop(columns=["sample_id", "phenotype", "ancestry"])
                
            if args.binary:
                res = run_gwas_logistic(X.to_numpy(), y, list(X.columns))
            else:
                res = run_gwas_linear(X.to_numpy(), y, list(X.columns))
    
            outfile = out_path.with_name(f"{out_path.stem}_{group}.csv")
            res.to_csv(outfile, index=False)
    
            print(f"Wrote ancestry-specific results: {outfile}")
