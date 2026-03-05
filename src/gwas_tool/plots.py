import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def manhattan_plot(df: pd.DataFrame, out_png, pval_col="p_value"):
    d = df.copy()
    d = d.dropna(subset=[pval_col])
    # Toy plot: x = SNP index (since toy dataset has no chr/pos yet)
    x = np.arange(d.shape[0])
    y = -np.log10(np.clip(d[pval_col].to_numpy(), 1e-300, 1.0))

    plt.figure()
    plt.scatter(x, y, s=10)
    plt.xlabel("SNP index")
    plt.ylabel("-log10(p)")
    plt.title("Manhattan (toy: SNP index)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def qq_plot(df: pd.DataFrame, out_png, pval_col="p_value"):
    d = df.copy()
    d = d.dropna(subset=[pval_col])
    p = np.sort(np.clip(d[pval_col].to_numpy(), 1e-300, 1.0))
    n = p.size
    exp = -np.log10(np.arange(1, n + 1) / (n + 1))
    obs = -np.log10(p)

    plt.figure()
    plt.scatter(exp, obs, s=10)
    maxv = max(exp.max(), obs.max())
    plt.plot([0, maxv], [0, maxv])
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title("QQ plot")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
