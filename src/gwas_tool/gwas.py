import numpy as np
from scipy import stats
import statsmodels.api as sm

def run_gwas_math(X, y, PCs=None):
    n, m = X.shape
    
    # Check if the phenotype is binary (0 and 1)
    is_binary = np.array_equal(y, y.astype(bool)) or set(np.unique(y)) <= {0, 1}

    if not is_binary:
        # --- LINEAR GWAS (Your original fast vectorized code) ---
        if PCs is not None:
            Z = np.hstack([np.ones((n, 1)), PCs])
            M_Z = np.eye(n) - Z @ np.linalg.inv(Z.T @ Z) @ Z.T
            y_star = M_Z @ y
            X_star = M_Z @ X
            df = n - Z.shape[1] - 1
        else:
            y_star = y - np.mean(y)
            X_star = X - np.mean(X, axis=0)
            df = n - 2
        
        numerator = np.dot(X_star.T, y_star)
        denominator = np.sum(X_star**2, axis=0)
        denominator[denominator == 0] = np.nan
        betas = numerator / denominator
        y_pred = X_star * betas
        rss = np.sum((y_star[:, np.newaxis] - y_pred)**2, axis=0)
        se = np.sqrt(rss / df) / np.sqrt(denominator)
        p_vals = stats.t.sf(np.abs(betas/se), df) * 2
        return betas, p_vals

    else:
        # --- LOGISTIC GWAS (New logic for binary traits) ---
        print("[*] Detecting binary phenotype. Running Logistic Regression...")
        betas = []
        p_vals = []
        
        # Prepare background covariates
        if PCs is not None:
            Z = sm.add_constant(PCs)
        else:
            Z = np.ones((n, 1))

        for i in range(m):
            snp_column = X[:, i]
            # Combine Intercept + PCs + SNP
            X_full = np.column_stack((Z, snp_column))
            
            try:
                model = sm.Logit(y, X_full)
                result = model.fit(disp=0)
                betas.append(result.params[-1])
                p_vals.append(result.pvalues[-1])
            except:
                betas.append(np.nan)
                p_vals.append(np.nan)
        
        return np.array(betas), np.array(p_vals)
