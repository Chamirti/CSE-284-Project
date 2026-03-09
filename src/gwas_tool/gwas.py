import pandas as pd
import numpy as np
import statsmodels.api as sm


def run_gwas(genotype_prefix, phenotype_file, covariate_file, output_file, model_type="linear"):
    """
    Run GWAS for each SNP using either linear or logistic regression.
    
    model_type:
        "linear"   -> continuous phenotype
        "logistic" -> binary phenotype (0/1)
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

        # Combine SNP + covariates
        X = np.hstack([snp_values, covariates])
        X = sm.add_constant(X)

        try:
            if model_type == "linear":
                model = sm.OLS(y, X).fit()

            elif model_type == "logistic":
                model = sm.Logit(y, X).fit(disp=0)

            else:
                raise ValueError("model_type must be 'linear' or 'logistic'")

            beta = model.params[1]
            pval = model.pvalues[1]

        except Exception:
            beta = np.nan
            pval = np.nan

        results.append({
            "SNP": snp,
            "beta": beta,
            "pval": pval
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
