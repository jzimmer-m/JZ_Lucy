from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from .auc import assert_columns


def read_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def resolve_path(base_dir: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base_dir / path).resolve()


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        return pd.read_excel(path)
    raise ValueError(f"Unsupported table extension for {path}")


def read_cyto_long_table(path: Path) -> pd.DataFrame:
    frame = read_table(path)
    assert_columns(
        frame,
        {
            "dataset",
            "context",
            "lnp",
            "time_h",
            "rep",
            "embedded_barcode",
            "cytosol_reads",
            "total_reads",
            "qpcr_norm_gl",
        },
        path,
    )
    return frame


def read_protein_long_table(
    path: Path,
    x_col: str = "time_h",
    value_col: str = "protein_expression",
) -> pd.DataFrame:
    frame = read_table(path)
    assert_columns(
        frame,
        {"dataset", "context", "lnp", x_col, "rep", value_col},
        path,
    )
    return frame


def read_color_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}

    frame = read_table(path)
    assert_columns(frame, {"lnp", "color"}, path)
    frame["lnp"] = frame["lnp"].astype(str).str.strip()
    frame["color"] = frame["color"].astype(str).str.strip()
    return dict(zip(frame["lnp"], frame["color"]))


def format_plex(value: object) -> str:
    if pd.isna(value):
        return ""

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value).strip()

    if numeric.is_integer():
        return str(int(numeric))
    return str(value).strip()


def prepare_lnp_metadata(
    rna: pd.DataFrame,
    keep_lnp_plexes: dict[str, list[object]] | None = None,
) -> pd.DataFrame:
    data = rna.copy()
    data["lnp"] = data["lnp"].astype(str).str.strip()
    data["include_lnp"] = True

    if "plex" not in data.columns:
        if keep_lnp_plexes:
            raise ValueError("RNA metadata needs a plex column to filter repeated LNPs")
        return data

    data["_plex_label"] = data["plex"].map(format_plex)

    if keep_lnp_plexes:
        for lnp, plexes in keep_lnp_plexes.items():
            lnp = str(lnp).strip()
            allowed = {format_plex(plex) for plex in plexes}
            is_excluded = data["lnp"].eq(lnp) & ~data["_plex_label"].isin(allowed)
            data.loc[is_excluded, "include_lnp"] = False

    unique_plexes = data[data["_plex_label"].ne("")].groupby("lnp")["_plex_label"].nunique()
    reused_lnps = unique_plexes[unique_plexes > 1].index
    is_reused = data["lnp"].isin(reused_lnps) & data["_plex_label"].ne("")
    data.loc[is_reused, "lnp"] = data.loc[is_reused, "lnp"] + " " + data.loc[is_reused, "_plex_label"]

    return data.drop(columns="_plex_label")


def apply_lnp_aliases(frame: pd.DataFrame, aliases: dict[str, str] | None) -> pd.DataFrame:
    if frame.empty or not aliases:
        return frame

    data = frame.copy()
    alias_map = {str(key).strip(): str(value).strip() for key, value in aliases.items()}
    data["lnp"] = data["lnp"].astype(str).str.strip().replace(alias_map)
    return data


def build_cyto_input(
    scar_df_path: Path,
    rna_metadata_path: Path,
    sample_metadata_path: Path,
    qpcr_path: Path,
    dataset: str,
    context_col: str,
    sample_join_cols: list[str],
    sample_filter: dict[str, str] | None = None,
    keep_lnp_plexes: dict[str, list[object]] | None = None,
) -> pd.DataFrame:
    scar = read_table(scar_df_path).rename(
        columns={
            "Cytosol": "cytosol_reads",
            "Endosome": "endosome_reads",
            "total": "total_reads",
        }
    )
    assert_columns(
        scar,
        {"flowcell", "barcode", "embedded_barcode", "cytosol_reads", "total_reads"},
        scar_df_path,
    )

    rna = read_table(rna_metadata_path)
    assert_columns(rna, {"internal_barcode", "LNP"}, rna_metadata_path)
    rna = rna.rename(columns={"internal_barcode": "embedded_barcode", "LNP": "lnp"})
    rna = prepare_lnp_metadata(rna, keep_lnp_plexes)

    samples = read_table(sample_metadata_path)
    assert_columns(samples, {"flowcell", "barcode", context_col, "time", "rep"}, sample_metadata_path)
    if sample_filter:
        for key, value in sample_filter.items():
            samples = samples[samples[key].astype(str).eq(str(value))]

    qpcr = read_table(qpcr_path)
    assert_columns(qpcr, set(sample_join_cols) | {"norm_GL"}, qpcr_path)

    joined_samples = samples.merge(qpcr, on=sample_join_cols, how="left", validate="many_to_one")
    missing_qpcr = joined_samples["norm_GL"].isna().sum()
    if missing_qpcr:
        raise ValueError(f"{missing_qpcr} samples in {sample_metadata_path} did not match qPCR rows")

    merged = (
        scar.merge(rna[["embedded_barcode", "lnp", "include_lnp"]], on="embedded_barcode", how="inner")
        .merge(joined_samples, on=["flowcell", "barcode"], how="inner")
        .rename(columns={context_col: "context", "time": "time_h", "norm_GL": "qpcr_norm_gl"})
    )
    merged["dataset"] = dataset
    return merged[
        [
            "dataset",
            "context",
            "lnp",
            "time_h",
            "rep",
            "embedded_barcode",
            "cytosol_reads",
            "total_reads",
            "qpcr_norm_gl",
            "include_lnp",
        ]
    ]
