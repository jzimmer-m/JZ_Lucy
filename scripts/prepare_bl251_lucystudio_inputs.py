from __future__ import annotations

import argparse
from pathlib import Path

from lucy_bl251.lucy_studio import normalize_lucy_studio_imaging_export


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize BL-251 LucyStudio protein-expression exports for AUC plotting.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/raw/LUCY-IVB-TEST-001.json"),
        help="LucyStudio project JSON exported from the 15 LNP Pilot.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/bl251_invitro_protein_expression.csv"),
        help="Tidy in vitro protein-expression CSV.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    normalize_lucy_studio_imaging_export(args.input.resolve()).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
