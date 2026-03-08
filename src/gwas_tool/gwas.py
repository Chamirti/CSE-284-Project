import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

def run_gwas_linear(genotype_prefix, phenotype_file, covariate_file, output_file):
    """
    Simple linear regression GWAS with covariates (PCs).
    Assumes PLINK .raw file format generated from genotype_prefix.
    """
    # Load genotype (PLINK .raw) and phenotype
    geno_file = genotype_prefix + ".raw"
    geno_df = pd.read_csv(geno_file, sep="\s+")
    pheno_df = pd.read_csv(phenotype_file, delim_whitespace=True, header=None)
    pcs_df = pd.read_csv(covariate_file, delim_whitespace=True, header=None)

    # Merge all
    data = pd.concat([geno_df, pheno_df, pcs_df], axis=1)
    y = data.iloc[:, -pheno_df.shape[1]:]  # phenotype column(s)
    X_cov = pcs_df.values                 # PCs as covariates

    results = []
    for snp in geno_df.columns[6:]:  # skip FID, IID, PAT, MAT, SEX, PHENOTYPE columns
        X_snp = geno_df[[snp]].values
        X = np.hstack([X_snp, X_cov])
        model = LinearRegression().fit(X, y)
        beta = model.coef_[0][0]
        residuals = y.values - model.predict(X)
        mse = np.mean(residuals**2)
        results.append({"SNP": snp, "Beta": beta, "MSE": mse})

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
