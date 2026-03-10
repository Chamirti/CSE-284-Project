import numpy as np
from sklearn.decomposition import PCA

def compute_pcs(geno_prefix, num_pcs=5):
    """
    Dummy PCA computation on genotype matrix
    """
    num_indivs = 50  # from your phenotype files
    num_snps = 1000
    np.random.seed(0)
    X = np.random.randint(0, 3, size=(num_indivs, num_snps))
    pca = PCA(n_components=num_pcs)
    pcs = pca.fit_transform(X)
    return pcs
