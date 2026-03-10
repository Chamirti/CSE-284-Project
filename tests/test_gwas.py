from src.gwas_tool.gwas import run_gwas_linear, run_gwas_logistic
from src.gwas_tool.pca import compute_pcs
from pathlib import Path

def test_linear():
    pcs = compute_pcs(None)
    pheno_file = Path("data/chr22/chr22_quantitative_clean.phen")
    results = run_gwas_linear(None, pheno_file, pcs)
    assert not results.empty

def test_logistic():
    pcs = compute_pcs(None)
    pheno_file = Path("data/chr22/chr22_casecontrol_clean.phen")
    results = run_gwas_logistic(None, pheno_file, pcs)
    assert not results.empty
