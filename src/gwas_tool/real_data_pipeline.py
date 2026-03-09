from pathlib import Path
import subprocess

from gwas_tool.gwas import run_gwas_linear
from gwas_tool.pca import compute_pcs
from gwas_tool.plots import manhattan_plot, qq_plot

def run_pipeline(geno_prefix, output_dir, gcta_path):
    geno_prefix = Path(geno_prefix)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Converting PLINK to RAW...")
    subprocess.run([
        "plink",
        "--bfile", str(geno_prefix),
        "--recode", "A",
        "--out", str(geno_prefix)
    ], check=True)

    print("Simulating phenotype with GCTA...")
    subprocess.run([
        gcta_path,
        "--bfile", str(geno_prefix),
        "--simu-qt", "1",
        "--out", str(output_dir / "simulated_pheno")
    ], check=True)

    pheno_file = output_dir / "simulated_pheno.phen"

    print("Computing PCA manually...")
    pcs_file = output_dir / "pcs.txt"
    compute_pcs(str(geno_prefix), pcs_file, num_pcs=10)

    print("Running GWAS...")
    gwas_results_file = output_dir / "gwas_results.csv"
    run_gwas_linear(
        genotype_prefix=str(geno_prefix),
        phenotype_file=str(pheno_file),
        covariate_file=str(pcs_file),
        output_file=str(gwas_results_file)
    )

    print("Generating plots...")
    bim_file = str(geno_prefix) + ".bim"
    manhattan_plot(str(gwas_results_file), bim_file, str(output_dir / "manhattan.png"))
    qq_plot(str(gwas_results_file), str(output_dir / "qq.png"))

    print("Pipeline complete!")
