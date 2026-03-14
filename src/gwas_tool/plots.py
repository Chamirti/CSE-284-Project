import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def generate_visuals(results, bim_path, causal_list, output_dir):
    bim = pd.read_csv(bim_path, sep=r'\s+', header=None, names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
    df = pd.merge(results, bim[['SNP', 'POS']], on='SNP')
    df['-log10P'] = -np.log10(df['P'])
    df = df.sort_values('POS')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

    ax1.scatter(df['POS'], df['-log10P'], c='royalblue', s=10, alpha=0.5, label='All SNPs')
    causal_df = df[df['SNP'].isin(causal_list)]
    ax1.scatter(causal_df['POS'], causal_df['-log10P'], c='red', s=40, label='Causal SNPs', edgecolors='black')
    ax1.axhline(y=-np.log10(5e-8), color='darkred', linestyle='--')
    ax1.set_title('Manhattan Plot')

    observed = np.sort(df['P'])
    expected = np.arange(1, len(observed) + 1) / (len(observed) + 1)
    ax2.scatter(-np.log10(expected), -np.log10(observed), c='black', s=10)
    max_val = np.max(-np.log10(expected))
    ax2.plot([0, max_val], [0, max_val], color='red', linestyle='--')
    ax2.set_title('QQ Plot')

    plt.tight_layout()
    plt.savefig(f"{output_dir}/gwas_visualizations.png")
