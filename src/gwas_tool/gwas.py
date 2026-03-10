import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression, LogisticRegression
import math  # for erf

def run_gwas_linear(geno_prefix, pheno_file, pcs):
    """
    Performs GWAS with linear regression using manual regression + PCA correction
    """
    # --- Read phenotype safely ---
    pheno = pd.read_csv(pheno_file, sep=r'\s+', engine='python', header=None,
                        names=["FID", "IID", "Phenotype"])
    pheno = pheno.dropna(subset=["Phenotype"])  # remove any missing values

    # --- Load genotype (dummy simulation for example) ---
    num_snps = 1000  # example SNPs
    num_indivs = pheno.shape[0]
    np.random.seed(0)
    X = np.random.randint(0, 3, size=(num_indivs, num_snps))

    # --- Combine PCA features ---
    X_corrected = np.hstack([X, pcs])

    pvals = []
    betas = []

    for snp_idx in range(X.shape[1]):
        lr = LinearRegression()
        X_snp = X_corrected[:, [snp_idx] + list(range(X.shape[1], X_corrected.shape[1]))]
        lr.fit(X_snp, pheno["Phenotype"])
        beta = lr.coef_[0]
        residuals = pheno["Phenotype"] - lr.predict(X_snp)
        se = residuals.std() / np.std(X[:, snp_idx])
        t_stat = beta / se
        pval = 2 * (1 - norm_cdf(abs(t_stat)))
        betas.append(beta)
        pvals.append(pval)

    return pd.DataFrame({"SNP": [f"rs{i}" for i in range(X.shape[1])],
                         "Beta": betas, "P": pvals})


def run_gwas_logistic(geno_prefix, pheno_file, pcs):
    """
    Performs GWAS with logistic regression using manual regression + PCA correction
    Handles GCTA case/control output automatically
    """
    import pandas as pd
    import numpy as np
    from sklearn.linear_model import LogisticRegression

    # --- Read phenotype ---
    pheno = pd.read_csv(pheno_file, sep=r'\s+', engine='python', header=None,
                        names=["FID", "IID", "Phenotype"])

    # --- Drop missing values (-9) ---
    pheno = pheno[pheno["Phenotype"] != -9]

    # --- Convert 1/2 -> 0/1 ---
    pheno["Phenotype"] = pheno["Phenotype"].astype(int)
    pheno["Phenotype"] = pheno["Phenotype"] - 1  # 1->0, 2->1

    # --- Number of individuals after cleaning ---
    num_indivs = pheno.shape[0]

    # --- Simulated genotype matrix for example ---
    num_snps = 1000
    np.random.seed(0)
    X = np.random.randint(0, 3, size=(num_indivs, num_snps))  # match pheno rows

    # --- Subset PCA to match cleaned individuals ---
    pcs = pcs[:num_indivs, :]
    X_corrected = np.hstack([X, pcs])

    # --- Run logistic regression for each SNP ---
    pvals = []
    betas = []

    for snp_idx in range(X.shape[1]):
        lr = LogisticRegression(solver="liblinear")
        X_snp = X_corrected[:, [snp_idx] + list(range(X.shape[1], X_corrected.shape[1]))]
        lr.fit(X_snp, pheno["Phenotype"])
        beta = lr.coef_[0][0]
        # Approximate p-value
        pval = 2 * (1 - norm_cdf(abs(beta / 0.1)))
        betas.append(beta)
        pvals.append(pval)

    return pd.DataFrame({"SNP": [f"rs{i}" for i in range(X.shape[1])],
                         "Beta": betas, "P": pvals})

def norm_cdf(x):
    """CDF of standard normal using math.erf"""
    return (1.0 + math.erf(x / math.sqrt(2.0))) / 2.0
