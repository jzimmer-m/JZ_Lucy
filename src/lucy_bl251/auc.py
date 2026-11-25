from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


CYTO_REQUIRED = {
    "dataset",
    "context",
    "lnp",
    "time_h",
    "rep",
    "embedded_barcode",
    "cytosol_reads",
    "total_reads",
    "qpcr_norm_gl",
}


def trapezoid_auc(x: pd.Series, y: pd.Series) -> float:
    """Return trapezoid AUC after sorting x/y pairs and dropping missing values."""
    points = (
        pd.DataFrame({"x": x, "y": y})
        .dropna()
        .sort_values("x")
        .drop_duplicates("x", keep="last")
    )
    if len(points) < 2:
        return math.nan
    return float(np.trapezoid(points["y"].to_numpy(), points["x"].to_numpy()))


def assert_columns(frame: pd.DataFrame, required: set[str], path: Path | str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def compute_cyto_abundance(frame: pd.DataFrame) -> pd.DataFrame:
    """Allocate qPCR-normalized GL abundance across LNP barcodes by cytosolic reads."""
    assert_columns(frame, CYTO_REQUIRED, "cyto_abund input")
    data = frame.copy()
    data["cytosol_reads"] = pd.to_numeric(data["cytosol_reads"], errors="coerce").fillna(0)
    data["total_reads"] = pd.to_numeric(data["total_reads"], errors="coerce").fillna(0)
    data["qpcr_norm_gl"] = pd.to_numeric(data["qpcr_norm_gl"], errors="coerce")

    sample_cols = ["dataset", "context", "time_h", "rep"]
    cytosol_totals = data.groupby(sample_cols)["cytosol_reads"].transform("sum")
    data["cyto_fraction"] = data["cytosol_reads"].div(cytosol_totals).where(cytosol_totals > 0, 0)
    data["cyto_abund"] = data["cyto_fraction"] * data["qpcr_norm_gl"]

    agg_spec = {
        "cyto_abund": ("cyto_abund", "sum"),
        "cytosol_reads": ("cytosol_reads", "sum"),
        "total_reads": ("total_reads", "sum"),
        "qpcr_norm_gl": ("qpcr_norm_gl", "first"),
    }
    if "include_lnp" in data.columns:
        agg_spec["include_lnp"] = ("include_lnp", "all")

    return data.groupby(["dataset", "context", "lnp", "time_h", "rep"], as_index=False).agg(
        **agg_spec
    )


def compute_auc(
    frame: pd.DataFrame,
    value_col: str,
    metric: str,
    dataset_col: str = "dataset",
    x_col: str = "time_h",
) -> pd.DataFrame:
    """Average replicates at each x value, then compute AUC for each dataset/context/LNP."""
    required = {dataset_col, "context", "lnp", x_col, value_col}
    assert_columns(frame, required, f"{metric} input")

    data = frame.copy()
    data[x_col] = pd.to_numeric(data[x_col], errors="coerce")
    data[value_col] = pd.to_numeric(data[value_col], errors="coerce")

    means = (
        data.dropna(subset=[x_col, value_col])
        .groupby([dataset_col, "context", "lnp", x_col], as_index=False)[value_col]
        .mean()
    )

    rows = []
    for keys, group in means.groupby([dataset_col, "context", "lnp"]):
        dataset, context, lnp = keys
        rows.append(
            {
                dataset_col: dataset,
                "context": context,
                "lnp": lnp,
                "auc": trapezoid_auc(group[x_col], group[value_col]),
            }
        )

    auc = pd.DataFrame(rows, columns=[dataset_col, "context", "lnp", "auc"])
    auc = auc.rename(columns={dataset_col: "dataset"})
    auc["metric"] = metric
    return auc[["metric", "dataset", "context", "lnp", "auc"]]


def make_comparison_matrix(auc: pd.DataFrame, value_col: str = "auc") -> pd.DataFrame:
    assert_columns(auc, {"metric", "dataset", "context", "lnp", value_col}, "comparison table")
    data = auc.copy()
    data["comparison"] = (
        data["metric"].astype(str)
        + ": "
        + data["dataset"].astype(str)
        + " - "
        + data["context"].astype(str)
    )
    return data.pivot_table(index="lnp", columns="comparison", values=value_col, aggfunc="mean")


def spearman_pair(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    pairs = pd.concat([x, y], axis=1).dropna()
    if len(pairs) < 2:
        return math.nan, math.nan, len(pairs)
    rho, p_value = stats.spearmanr(pairs.iloc[:, 0], pairs.iloc[:, 1])
    return float(rho), float(p_value), len(pairs)
