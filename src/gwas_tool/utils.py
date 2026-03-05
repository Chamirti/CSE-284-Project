from __future__ import annotations
import numpy as np
import pandas as pd


def load_geno_csv(path: str) -> pd.DataFrame:
    """
    Load genotype CSV with first column `sample_id` and remaining columns SNPs coded as 0/1/2.
    """
    df = pd.read_csv(path)
    if "sample_id" not in df.columns:
        raise ValueError("Genotype file must contain a `sample_id` column.")
    return df


def load_pheno_csv(path: str) -> pd.DataFrame:
    """
    Load phenotype CSV with columns: `sample_id`, `phenotype`.
    """
    df = pd.read_csv(path)
    required = {"sample_id", "phenotype"}
    if not required.issubset(set(df.columns)):
        raise ValueError("Phenotype file must contain columns: sample_id, phenotype.")
    return df[["sample_id", "phenotype"]]


def align_geno_pheno(geno: pd.DataFrame, pheno: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """
    Inner join geno and pheno on sample_id.
    Returns:
      X_df: genotype matrix (samples x SNPs) as DataFrame
      y: phenotype vector (n,)
      snp_names: list of SNP column names
    """
    merged = pheno.merge(geno, on="sample_id", how="inner")
    if merged.shape[0] == 0:
        raise ValueError("No overlapping sample_id between genotype and phenotype files.")

    y = merged["phenotype"].to_numpy()
    X_df = merged.drop(columns=["sample_id", "phenotype"])
    snp_names = list(X_df.columns)
    return X_df, y, snp_names


def clip_pvalues(p: np.ndarray, eps: float = 1e-300) -> np.ndarray:
    """
    Avoid log(0) in plots.
    """
    p = np.asarray(p, dtype=float)
    return np.clip(p, eps, 1.0)
