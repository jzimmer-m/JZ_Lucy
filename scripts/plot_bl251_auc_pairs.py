from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from lucy_bl251.auc import compute_auc, compute_cyto_abundance
from lucy_bl251.io import (
    apply_lnp_aliases,
    build_cyto_input,
    read_color_map,
    read_config,
    read_cyto_long_table,
    read_protein_long_table,
    resolve_path,
)
from lucy_bl251.plot import plot_auc_pairs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare BL-251 Scar-seq time AUCs with protein expression at 40 ng.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/bl251_auc_inputs.yml"),
        help="YAML file describing Scar-seq, protein, and color inputs.",
    )
    return parser.parse_args()


def load_cyto_inputs(config: dict, base_dir: Path) -> pd.DataFrame:
    frames = []
    for item in config.get("scarseq", []):
        if "cyto_abund" in item:
            frame = read_cyto_long_table(resolve_path(base_dir, item["cyto_abund"]))
        else:
            frame = build_cyto_input(
                scar_df_path=resolve_path(base_dir, item["scar_df"]),
                rna_metadata_path=resolve_path(base_dir, item["rna_metadata"]),
                sample_metadata_path=resolve_path(base_dir, item["sample_metadata"]),
                qpcr_path=resolve_path(base_dir, item["qpcr"]),
                dataset=item["dataset"],
                context_col=item["context_col"],
                sample_join_cols=item["sample_join_cols"],
                sample_filter=item.get("sample_filter"),
                keep_lnp_plexes=item.get("keep_lnp_plexes"),
            )
        cyto_abund = compute_cyto_abundance(frame)
        if "include_lnp" in cyto_abund:
            cyto_abund = cyto_abund[cyto_abund["include_lnp"]].drop(columns="include_lnp")
        frames.append(cyto_abund)
    return pd.concat(frames, ignore_index=True)


def load_protein_40ng_values(config: dict, base_dir: Path) -> list[pd.DataFrame]:
    value_frames = []
    for item in config.get("protein", []):
        dose = float(item.get("dose_ng", 40))
        value_col = item.get("value_col", "protein_expression")
        metric = item.get("metric", f"In vitro protein expression at {dose:g} ng")
        frame = read_protein_long_table(
            resolve_path(base_dir, item["protein_expression"]),
            x_col="dose_ng",
            value_col=value_col,
        )
        frame = apply_lnp_aliases(frame, config.get("lnp_aliases"))
        frame["dose_ng"] = pd.to_numeric(frame["dose_ng"], errors="coerce")
        frame[value_col] = pd.to_numeric(frame[value_col], errors="coerce")
        selected = frame.loc[frame["dose_ng"].eq(dose)].dropna(subset=[value_col])
        if selected.empty:
            raise ValueError(f"No protein-expression measurements at {dose:g} ng")
        values = (
            selected.groupby(["dataset", "context", "lnp"], as_index=False)[value_col]
            .mean()
            .rename(columns={value_col: "value"})
        )
        values["metric"] = metric
        values["measure"] = "single_dose_mean"
        value_frames.append(values[["metric", "dataset", "context", "lnp", "measure", "value"]])
    return value_frames


def main() -> None:
    args = parse_args()
    repo_dir = args.config.resolve().parent.parent
    config = read_config(args.config)

    cyto = apply_lnp_aliases(load_cyto_inputs(config, repo_dir), config.get("lnp_aliases"))
    scar_auc = compute_auc(cyto, "cyto_abund", "Scar-seq cyto abundance")
    auc_path = resolve_path(repo_dir, config["outputs"]["auc_csv"])
    auc_path.parent.mkdir(parents=True, exist_ok=True)
    scar_auc.to_csv(auc_path, index=False)

    scar_values = scar_auc.rename(columns={"auc": "value"}).copy()
    scar_values["measure"] = "time_auc"
    values = pd.concat([scar_values, *load_protein_40ng_values(config, repo_dir)], ignore_index=True)
    values = values[["metric", "dataset", "context", "lnp", "measure", "value"]]
    comparison_path = resolve_path(repo_dir, config["outputs"]["comparison_csv"])
    comparison_path.parent.mkdir(parents=True, exist_ok=True)
    values.to_csv(comparison_path, index=False)

    color_path = config.get("lnp_colors")
    color_map = read_color_map(resolve_path(repo_dir, color_path) if color_path else None)
    figure_defaults = config.get("figure_style", {})
    figures = config["outputs"].get(
        "figures",
        [
            {
                "figure_pdf": config["outputs"]["figure_pdf"],
                "figure_png": config["outputs"].get("figure_png"),
                "comparison_pairs": config.get("comparison_pairs"),
            }
        ],
    )
    for figure in figures:
        figure_style = {**figure_defaults, **figure.get("style", {})}
        plot_auc_pairs(
            auc=values,
            output_pdf=resolve_path(repo_dir, figure["figure_pdf"]),
            output_png=resolve_path(repo_dir, figure.get("figure_png"))
            if figure.get("figure_png")
            else None,
            color_map=color_map,
            comparison_pairs=figure.get("comparison_pairs", config.get("comparison_pairs")),
            **figure_style,
        )


if __name__ == "__main__":
    main()
