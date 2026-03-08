import subprocess
from pathlib import Path
from gwas import run_gwas_linear
from pca import compute_pcs
from plots import manhattan_plot, qq_plot

# ------------------------------
# File paths
# ------------------------------
geno_prefix = Path("data/chr22/chr22")
gcta_path = Path("/path/to/gcta")
output_dir = Path("results/real_data")
output_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------
# Step 1: Simulate phenotype
# ------------------------------
pheno_file = output_dir / "simulated_pheno.txt"
subprocess.run([
    str(gcta_path),
    "--bfile", str(geno_prefix),
    "--simu-qt", "1",
    "--out", str(output_dir / "simulated_pheno")
])
print(f"Phenotypes saved to {pheno_file}")

# ------------------------------
# Step 2: PCA
# ------------------------------
pcs_file = output_dir / "pcs.txt"
compute_pcs(str(geno_prefix), str(pcs_file))
print(f"PCA saved to {pcs_file}")

# ------------------------------
# Step 3: GWAS
# ------------------------------
gwas_results_file = output_dir / "gwas_results.csv"
run_gwas_linear(
    genotype_prefix=str(geno_prefix),
    phenotype_file=str(pheno_file),
    covariate_file=str(pcs_file),
    output_file=str(gwas_results_file)
)
print(f"GWAS results saved to {gwas_results_file}")

# ------------------------------
# Step 4: Plots
# ------------------------------
manhattan_plot(str(gwas_results_file), str(output_dir / "manhattan.png"))
qq_plot(str(gwas_results_file), str(output_dir / "qq.png"))
print("Pipeline complete!")
