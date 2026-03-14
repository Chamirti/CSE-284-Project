import numpy as np
from scipy import stats

def run_gwas_math(X, y, PCs=None):
    n, m = X.shape
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
