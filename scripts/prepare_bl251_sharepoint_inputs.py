from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from lucy_bl251.sharepoint_results import normalize_scarseq_result_export


CYTO_EXPORTS = [
    {
        "dataset": "in vivo Scar-seq",
        "input_name": "LUCY-IVM-Results-BL251Scarseq-2026-11-20-JZ.csv",
        "output_name": "bl251_invivo_cyto_abundance.csv",
        "context_col": "tissue",
        "exclude_contexts": (),
    },
    {
        "dataset": "in vitro Scar-seq",
        "input_name": "LUCY-IVM-Results-invitroBL251Scarseq-2026-11-20-JZ.csv",
        "output_name": "bl251_invitro_cyto_abundance.csv",
        "context_col": "cell",
        "exclude_contexts": ("DS",),
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize BL-251 Project Lucy SharePoint CSV exports for AUC plotting.",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("data/raw"),
        help="Directory containing the raw SharePoint result CSVs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed"),
        help="Directory for normalized cyto-abundance tables.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw_dir = args.raw_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for export in CYTO_EXPORTS:
        frame = pd.read_csv(raw_dir / export["input_name"])
        normalized = normalize_scarseq_result_export(
            frame,
            dataset=export["dataset"],
            context_col=export["context_col"],
            exclude_contexts=export["exclude_contexts"],
            keep_lnp_plexes={"86/56": [20]},
        )
        normalized.to_csv(output_dir / export["output_name"], index=False)


if __name__ == "__main__":
    main()
