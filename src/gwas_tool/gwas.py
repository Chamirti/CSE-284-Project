import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression, LogisticRegression

def run_gwas_linear(geno_prefix, pheno_file, pcs):
    """
    Performs GWAS with linear regression using manual regression + PCA correction
    """
    # Load phenotype
    pheno = pd.read_csv(pheno_file, sep=" ", header=None, names=["FID","IID","Phenotype"])
    
    # Load genotype (dummy simulation for example)
    num_snps = 1000  # example number
    num_indivs = pheno.shape[0]
    np.random.seed(0)
    X = np.random.randint(0, 3, size=(num_indivs, num_snps))

    # Combine PCA features
    X_corrected = np.hstack([X, pcs])

    pvals = []
    betas = []

    for snp_idx in range(X.shape[1]):
        lr = LinearRegression()
        lr.fit(X_corrected[:, [snp_idx]+list(range(X.shape[1], X_corrected.shape[1]))], pheno["Phenotype"])
        beta = lr.coef_[0]
        residuals = pheno["Phenotype"] - lr.predict(X_corrected[:, [snp_idx]+list(range(X.shape[1], X_corrected.shape[1]))])
        se = residuals.std() / np.std(X[:, snp_idx])
        t_stat = beta / se
        pval = 2*(1 - norm_cdf(np.abs(t_stat)))
        betas.append(beta)
        pvals.append(pval)

    return pd.DataFrame({"SNP": [f"rs{i}" for i in range(X.shape[1])], "Beta": betas, "P": pvals})

def run_gwas_logistic(geno_prefix, pheno_file, pcs):
    """
    Performs GWAS with logistic regression using manual regression + PCA correction
    """
    pheno = pd.read_csv(pheno_file, sep=" ", header=None, names=["FID","IID","Phenotype"])
    
    num_snps = 1000
    num_indivs = pheno.shape[0]
    np.random.seed(0)
    X = np.random.randint(0, 3, size=(num_indivs, num_snps))
    
    X_corrected = np.hstack([X, pcs])

    pvals = []
    betas = []

    for snp_idx in range(X.shape[1]):
        lr = LogisticRegression(solver="liblinear")
        lr.fit(X_corrected[:, [snp_idx]+list(range(X.shape[1], X_corrected.shape[1]))], pheno["Phenotype"])
        beta = lr.coef_[0][0]
        # For simplicity, approximate p-value
        pval = 2*(1 - norm_cdf(np.abs(beta/0.1)))
        betas.append(beta)
        pvals.append(pval)

    return pd.DataFrame({"SNP": [f"rs{i}" for i in range(X.shape[1])], "Beta": betas, "P": pvals})

def norm_cdf(x):
    return (1.0 + np.math.erf(x / np.sqrt(2.0))) / 2.0
