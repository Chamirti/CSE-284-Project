import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression


def _add_intercept_and_cov(X_snp: np.ndarray, covariates: np.ndarray | None):
    # X_snp: (n,) single SNP
    n = X_snp.shape[0]
    intercept = np.ones((n, 1))
    snp_col = X_snp.reshape(n, 1)

    if covariates is None:
        return np.hstack([intercept, snp_col])
    return np.hstack([intercept, snp_col, covariates])


def run_gwas_linear(G: np.ndarray, y: np.ndarray, snp_names: list[str], covariates: np.ndarray | None = None) -> pd.DataFrame:
    """
    Linear regression per SNP:
      y = b0 + b1*SNP + cov + e

    Returns beta and p-value for SNP coefficient (b1).
    """
    y = y.astype(float)
    n, m = G.shape
    results = []

    for j in range(m):
        X = _add_intercept_and_cov(G[:, j], covariates)
        # OLS closed form: b = (X'X)^-1 X'y
        XtX = X.T @ X
        try:
            XtX_inv = np.linalg.inv(XtX)
        except np.linalg.LinAlgError:
            # singular; skip
            results.append((snp_names[j], np.nan, np.nan))
            continue
        b = XtX_inv @ (X.T @ y)
        y_hat = X @ b
        resid = y - y_hat

        dof = X.shape[0] - X.shape[1]
        if dof <= 0:
            results.append((snp_names[j], float(b[1]), np.nan))
            continue

        s2 = (resid @ resid) / dof
        se_b = np.sqrt(np.diag(XtX_inv) * s2)

        t_stat = b[1] / se_b[1] if se_b[1] > 0 else np.nan
        p_val = 2 * stats.t.sf(np.abs(t_stat), df=dof) if np.isfinite(t_stat) else np.nan

        results.append((snp_names[j], float(b[1]), float(p_val)))

    return pd.DataFrame(results, columns=["snp", "beta", "p_value"])


def run_gwas_logistic(G: np.ndarray, y: np.ndarray, snp_names: list[str], covariates: np.ndarray | None = None) -> pd.DataFrame:
    """
    Logistic regression per SNP using sklearn LogisticRegression (L2, large C ~ unregularized).
    p-values approximated via Wald z-test from fitted coefficients and (approx) Hessian.

    Note: For the class project (toy data), this is acceptable; later you can replace with
    statsmodels for exact SE/p-values if desired.
    """
    y = y.astype(int)
    if set(np.unique(y)) - {0, 1}:
        raise ValueError("Binary phenotype must be 0/1 for logistic GWAS.")

    n, m = G.shape
    results = []

    for j in range(m):
        X = _add_intercept_and_cov(G[:, j], covariates)

        # Fit logistic regression; make C very large to reduce regularization effect
        model = LogisticRegression(
            penalty="l2",
            C=1e6,
            solver="lbfgs",
            max_iter=1000
        )
        try:
            model.fit(X[:, 1:], y)  # exclude intercept because sklearn adds its own intercept by default if fit_intercept=True
        except Exception:
            results.append((snp_names[j], np.nan, np.nan))
            continue

        # sklearn stores intercept separately; coefficient for SNP is coef_[0,0] because our first column is SNP when intercept removed
        beta_snp = float(model.coef_[0][0])

        # Approximate standard error using observed information matrix:
        # For logistic regression: I = X^T W X, where W = p*(1-p)
        # Use X_design without intercept (since model handled intercept), approximate around fitted probs.
        p_hat = model.predict_proba(X[:, 1:])[:, 1]
        W = p_hat * (1 - p_hat)
        Xd = X[:, 1:]  # includes SNP + covariates (no intercept)
        XtWX = Xd.T @ (Xd * W.reshape(-1, 1))
        try:
            cov_mat = np.linalg.inv(XtWX)
            se = float(np.sqrt(cov_mat[0, 0]))  # SNP term
            z = beta_snp / se if se > 0 else np.nan
            p_val = 2 * stats.norm.sf(np.abs(z)) if np.isfinite(z) else np.nan
        except np.linalg.LinAlgError:
            p_val = np.nan

        results.append((snp_names[j], beta_snp, float(p_val)))

    return pd.DataFrame(results, columns=["snp", "beta", "p_value"])
