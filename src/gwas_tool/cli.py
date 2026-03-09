import argparse
from pathlib import Path
from gwas_tool.real_data_pipeline import run_pipeline


def main():

    parser = argparse.ArgumentParser(
        description="GWAS Tool: Run a GWAS pipeline on genotype data"
    )

    parser.add_argument(
        "--geno",
        required=True,
        help="Prefix of PLINK genotype files (.bed/.bim/.fam)"
    )

    parser.add_argument(
        "--out",
        default="results",
        help="Output directory"
    )

    parser.add_argument(
        "--gcta",
        default="gcta64",
        help="Path to GCTA executable"
    )

    args = parser.parse_args()

    run_pipeline(
        geno_prefix=args.geno,
        output_dir=args.out,
        gcta_path=args.gcta
    )


if __name__ == "__main__":
    main()
