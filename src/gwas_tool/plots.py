import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def manhattan_plot(gwas_results_file, output_file):
    df = pd.read_csv(gwas_results_file)
    plt.figure(figsize=(12,6))
    plt.scatter(range(len(df)), -np.log10(df['MSE']), c='blue', s=10)
    plt.xlabel("SNPs")
    plt.ylabel("-log10(MSE)")
    plt.title("Manhattan Plot")
    plt.savefig(output_file)
    plt.close()

def qq_plot(gwas_results_file, output_file):
    df = pd.read_csv(gwas_results_file)
    expected = -np.log10(np.linspace(1/len(df),1,len(df)))
    observed = -np.log10(np.sort(df['MSE']))
    plt.figure(figsize=(6,6))
    plt.scatter(expected, observed)
    plt.plot([0,max(expected)], [0,max(expected)], color='red', linestyle='--')
    plt.xlabel("Expected -log10(MSE)")
    plt.ylabel("Observed -log10(MSE)")
    plt.title("QQ Plot")
    plt.savefig(output_file)
    plt.close()
