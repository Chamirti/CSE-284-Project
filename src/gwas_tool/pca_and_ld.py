import numpy as np

def run_ld_pruning(X, window_size=50, step_size=5, r2_threshold=0.2):
    n, m = X.shape
    kept_indices = []
    for i in range(0, m - window_size, step_size):
        window = X[:, i:i + window_size]
        corr_matrix = np.corrcoef(window, rowvar=False)
        for j in range(window_size):
            idx = i + j
            is_redundant = False
            for k in range(j):
                if corr_matrix[j, k]**2 > r2_threshold:
                    is_redundant = True
                    break
            if not is_redundant and idx not in kept_indices:
                kept_indices.append(idx)
    return kept_indices

def run_pca(X_pruned, n_pcs=3):
    X_std_pruned = (X_pruned - np.mean(X_pruned, axis=0)) / np.std(X_pruned, axis=0)
    X_std_pruned = np.nan_to_num(X_std_pruned)
    U, S, Vt = np.linalg.svd(X_std_pruned, full_matrices=False)
    return U[:, :n_pcs]
