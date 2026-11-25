from __future__ import annotations

from collections.abc import Collection, Mapping

import pandas as pd

from .auc import assert_columns
from .io import format_plex


CYTO_EXPORT_COLUMNS = {
    "internal_barcode",
    "Cytosol",
    "total",
    "LNP",
    "time",
    "rep",
    "norm_GL",
    "plex",
}

CYTO_LONG_COLUMNS = [
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


def normalize_scarseq_result_export(
    frame: pd.DataFrame,
    *,
    dataset: str,
    context_col: str,
    exclude_contexts: Collection[str] = (),
    keep_lnp_plexes: Mapping[str, Collection[object]] | None = None,
) -> pd.DataFrame:
    assert_columns(frame, CYTO_EXPORT_COLUMNS | {context_col}, "SharePoint Scar-seq export")

    data = frame.rename(
        columns={
            context_col: "context",
            "LNP": "lnp",
            "time": "time_h",
            "internal_barcode": "embedded_barcode",
            "Cytosol": "cytosol_reads",
            "total": "total_reads",
            "norm_GL": "qpcr_norm_gl",
        }
    ).copy()
    data["dataset"] = dataset
    data["lnp"] = data["lnp"].astype(str).str.strip()

    if exclude_contexts:
        excluded = {str(context) for context in exclude_contexts}
        data = data[~data["context"].astype(str).isin(excluded)].copy()

    data = mark_repeated_lnp_plexes(data, keep_lnp_plexes)
    data["time_h"] = pd.to_numeric(data["time_h"], errors="coerce")

    return data[CYTO_LONG_COLUMNS].reset_index(drop=True)


def mark_repeated_lnp_plexes(
    frame: pd.DataFrame,
    keep_lnp_plexes: Mapping[str, Collection[object]] | None = None,
) -> pd.DataFrame:
    data = frame.copy()
    data["include_lnp"] = True
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
