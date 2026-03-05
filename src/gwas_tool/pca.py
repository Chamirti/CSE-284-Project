import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def compute_pcs(geno_df: pd.DataFrame, n_components: int = 2) -> pd.DataFrame:
    """
    Compute PCs from genotype matrix (samples x SNPs).
    Standardize SNP columns before PCA.
    """
    X = geno_df.to_numpy()
    Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    pcs = PCA(n_components=n_components, random_state=0).fit_transform(Xs)
    cols = [f"PC{i+1}" for i in range(n_components)]
    return pd.DataFrame(pcs, columns=cols)
