import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def manhattan_plot(gwas_results_file, bim_file, output_file):
    gwas = pd.read_csv(gwas_results_file)
    bim = pd.read_csv(bim_file, delim_whitespace=True, header=None,
                      names=["CHR","SNP","CM","BP","A1","A2"])
    merged = gwas.merge(bim, on="SNP")
    merged["minus_log10"] = -np.log10(merged["pval"] + 1e-12)
    merged = merged.sort_values(["CHR","BP"])

    plt.figure(figsize=(12,6))
    colors = ["#4C72B0", "#55A868"]
    x_vals, tick_positions, tick_labels = [], [], []
    last_base = 0

    for i, (chrom, group) in enumerate(merged.groupby("CHR")):
        x = group["BP"] + last_base
        plt.scatter(x, group["minus_log10"], c=colors[i % 2], s=8)
        x_vals.extend(x)
        tick_positions.append(x.mean())
        tick_labels.append(str(chrom))
        last_base = x.max()

    sig_line = -np.log10(5e-8)
    plt.axhline(sig_line, color="red", linestyle="--")
    plt.xticks(tick_positions, tick_labels)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(p-value)")
    plt.title("Manhattan Plot")
    plt.savefig(output_file)
    plt.close()

def qq_plot(gwas_results_file, output_file):
    df = pd.read_csv(gwas_results_file)
    observed = -np.log10(np.sort(df["pval"] + 1e-12))
    expected = -np.log10(np.linspace(1/len(observed),1,len(observed)))
    plt.figure(figsize=(6,6))
    plt.scatter(expected, observed)
    plt.plot([0,max(expected)], [0,max(expected)], linestyle="--")
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title("QQ Plot")
    plt.savefig(output_file)
    plt.close()
