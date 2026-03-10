import argparse
from pathlib import Path
from .gwas import run_gwas_linear, run_gwas_logistic
from .pca import compute_pcs
from .plots import manhattan_plot, qq_plot

def main():
    parser = argparse.ArgumentParser(description="GWAS Tool: Run manual GWAS with PCA correction")
    parser.add_argument("--geno", type=str, required=True, help="Path prefix to genotype file (PLINK .bed/.bim/.fam)")
    parser.add_argument("--pheno", type=str, required=True, help="Path to phenotype file")
    parser.add_argument("--type", type=str, choices=["linear", "logistic"], required=True, help="Regression type")
    parser.add_argument("--out", type=str, required=True, help="Output directory")

    args = parser.parse_args()

    geno_prefix = Path(args.geno)
    pheno_file = Path(args.pheno)
    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True, parents=True)

    print("Computing PCA...")
    pcs = compute_pcs(geno_prefix)
    
    print(f"Running {args.type} GWAS...")
    if args.type == "linear":
        results = run_gwas_linear(geno_prefix, pheno_file, pcs)
    else:
        results = run_gwas_logistic(geno_prefix, pheno_file, pcs)

    results_file = out_dir / f"gwas_{args.type}_results.csv"
    results.to_csv(results_file, index=False)
    print(f"Results saved to {results_file}")

    print("Generating plots...")
    manhattan_plot(results, out_dir / f"manhattan_{args.type}.png")
    qq_plot(results, out_dir / f"qq_{args.type}.png")
    print("Done.")

if __name__ == "__main__":
    main()
