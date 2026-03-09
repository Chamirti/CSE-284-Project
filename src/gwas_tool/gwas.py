import pandas as pd
import numpy as np
import statsmodels.api as sm

def run_gwas_linear(genotype_prefix, phenotype_file, covariate_file, output_file):
    """
    Run linear GWAS for each SNP, correcting for covariates (PCs).
    """
    geno_file = genotype_prefix + ".raw"
    geno_df = pd.read_csv(geno_file, delim_whitespace=True)
    pheno_df = pd.read_csv(phenotype_file, delim_whitespace=True, header=None)
    pcs_df = pd.read_csv(covariate_file, delim_whitespace=True, header=None)

    y = pheno_df.iloc[:, -1].values
    covariates = pcs_df.values

    results = []

    for snp in geno_df.columns[6:]:
        snp_values = geno_df[snp].values.reshape(-1, 1)
        X = np.hstack([snp_values, covariates])
        X = sm.add_constant(X)
        model = sm.OLS(y, X).fit()
        beta = model.params[1]
        pval = model.pvalues[1]
        results.append({"SNP": snp, "beta": beta, "pval": pval})

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
