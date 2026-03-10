import matplotlib.pyplot as plt
import numpy as np

def manhattan_plot(results, out_file):
    plt.figure(figsize=(12,6))
    plt.scatter(range(len(results)), -np.log10(results["P"]), c='blue', s=10)
    plt.xlabel("SNP index")
    plt.ylabel("-log10(p-value)")
    plt.title("Manhattan Plot")
    plt.savefig(out_file)
    plt.close()

def qq_plot(results, out_file):
    pvals_sorted = np.sort(results["P"])
    expected = np.arange(1, len(pvals_sorted)+1)/len(pvals_sorted)
    plt.figure(figsize=(6,6))
    plt.scatter(-np.log10(expected), -np.log10(pvals_sorted))
    plt.plot([0, -np.log10(expected[0])], [0, -np.log10(expected[0])], color='red')
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.title("QQ Plot")
    plt.savefig(out_file)
    plt.close()
