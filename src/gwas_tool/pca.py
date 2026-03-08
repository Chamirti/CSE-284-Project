import subprocess

def compute_pcs(genotype_prefix, output_file, num_pcs=10):
    """
    Uses PLINK to compute top principal components for population structure.
    """
    # Run PLINK PCA
    subprocess.run([
        "plink",
        "--bfile", genotype_prefix,
        "--pca", str(num_pcs),
        "--out", output_file.with_suffix('')
    ])
