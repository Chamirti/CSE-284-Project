import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

def compute_pcs(genotype_prefix, output_file, num_pcs=10):
    """
    Compute PCA manually from PLINK RAW file.
    """
    raw_file = str(genotype_prefix) + ".raw"
    geno_df = pd.read_csv(raw_file, delim_whitespace=True)

    snp_data = geno_df.iloc[:, 6:].values

    # Standardize SNPs
    snp_mean = snp_data.mean(axis=0)
    snp_std = snp_data.std(axis=0, ddof=1)
    snp_std[snp_std == 0] = 1
    snp_scaled = (snp_data - snp_mean) / snp_std

    # PCA
    pca = PCA(n_components=num_pcs)
    pcs = pca.fit_transform(snp_scaled)

    pcs_df = pd.DataFrame(
        np.hstack([geno_df.iloc[:, :2].values, pcs]),
        columns=["FID", "IID"] + [f"PC{i+1}" for i in range(num_pcs)]
    )

    pcs_df.to_csv(output_file, sep=" ", index=False, header=False)
    print("Explained variance by PCs:", pca.explained_variance_ratio_)
