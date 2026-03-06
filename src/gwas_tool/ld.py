import numpy as np
import pandas as pd


def compute_r2(x: np.ndarray, y: np.ndarray) -> float:
    
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.size != y.size:
        raise ValueError("Genotype vectors must have the same length.")

    # If either SNP is constant, correlation is undefined
    if np.var(x) == 0 or np.var(y) == 0:
        return 0.0

    r = np.corrcoef(x, y)[0, 1]
    if np.isnan(r):
        return 0.0
    return float(r * r)


def get_neighbor_indices(lead_idx: int, n_snps: int, window_snps: int) -> list[int]:
   
    start = max(0, lead_idx - window_snps)
    end = min(n_snps, lead_idx + window_snps + 1)

    neighbors = []
    for idx in range(start, end):
        if idx != lead_idx:
            neighbors.append(idx)
    return neighbors


def ld_clump(
    geno_df: pd.DataFrame,
    gwas_df: pd.DataFrame,
    ld_threshold: float = 0.8,
    window_snps: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    
    required_cols = {"snp", "p_value"}
    if not required_cols.issubset(gwas_df.columns):
        raise ValueError("gwas_df must contain columns: snp, p_value")

    snp_names = list(geno_df.columns)
    snp_to_idx = {snp: i for i, snp in enumerate(snp_names)}
    n_snps = len(snp_names)

    # Keep only SNPs that appear in both genotype and GWAS results
    work_df = gwas_df[gwas_df["snp"].isin(snp_names)].copy()
    work_df = work_df.dropna(subset=["p_value"])
    work_df = work_df.sort_values("p_value", ascending=True).reset_index(drop=True)

    remaining = set(work_df["snp"].tolist())
    lead_rows = []
    clump_rows = []

    for _, row in work_df.iterrows():
        lead_snp = row["snp"]
        lead_p = row["p_value"]

        if lead_snp not in remaining:
            continue

        lead_idx = snp_to_idx[lead_snp]
        lead_vec = geno_df[lead_snp].to_numpy()

        # This SNP becomes a lead SNP
        clumped_members = []
        neighbor_indices = get_neighbor_indices(lead_idx, n_snps, window_snps)

        for nb_idx in neighbor_indices:
            neighbor_snp = snp_names[nb_idx]

            if neighbor_snp not in remaining:
                continue

            neighbor_vec = geno_df[neighbor_snp].to_numpy()
            r2 = compute_r2(lead_vec, neighbor_vec)

            if r2 >= ld_threshold:
                remaining.remove(neighbor_snp)
                clumped_members.append((neighbor_snp, r2))

       
        remaining.remove(lead_snp)

        lead_rows.append(
            {
                "lead_snp": lead_snp,
                "p_value": lead_p,
                "num_clumped": len(clumped_members),
            }
        )

        for member_snp, r2 in clumped_members:
            clump_rows.append(
                {
                    "lead_snp": lead_snp,
                    "member_snp": member_snp,
                    "r2": r2,
                }
            )

    lead_df = pd.DataFrame(lead_rows)
    clump_df = pd.DataFrame(clump_rows)

    return lead_df, clump_df
