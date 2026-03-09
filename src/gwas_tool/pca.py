import subprocess
import pandas as pd


def compute_pcs(genotype_prefix, output_file, num_pcs=10):

    prefix = str(output_file).replace(".txt","")

    subprocess.run([
        "plink",
        "--bfile", genotype_prefix,
        "--pca", str(num_pcs),
        "--out", prefix
    ])

    eigenvec_file = prefix + ".eigenvec"

    df = pd.read_csv(eigenvec_file, delim_whitespace=True, header=None)

    pcs = df.iloc[:,2:]   # remove FID IID

    pcs.to_csv(output_file, sep=" ", index=False, header=False)
