from __future__ import annotations

from collections.abc import Iterable, Mapping
from itertools import combinations, product
from math import ceil
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullLocator
import numpy as np
import pandas as pd
import seaborn as sns

from .auc import make_comparison_matrix, spearman_pair


DEFAULT_COLORS = {
    "86/56": "#4C78A8",
    "716/56": "#F58518",
    "2210/56": "#E45756",
    "2210/DMPS": "#EECA3B",
    "86/56 low C:D": "#FF9DA6",
    "86/56 +A4h4 (LNP-20A)": "#AF7AA1",
    "2210/DSG": "#00A6A6",
    "2210/DOPE/DSG": "#5B8FF9",
    "7628/PG": "#6C8E3F",
    "86/56 Legacy": "#8C564B",
    "86/56 PIPA 2 on RAIDD": "#D62728",
}


def comparison_label(comparison: str) -> str:
    suffix = "" if "protein expression at 40 ng" in comparison.lower() else " AUC"
    return fill(f"{comparison}{suffix}", width=36)


def axis_label(comparisons: list[str], metadata: Mapping[str, Mapping[str, str]]) -> str:
    first = metadata[comparisons[0]]
    if first["dataset"] == "in vitro protein expression":
        return "In vitro protein expression at 40 ng, 24 h (Alexa 488 mean per well)"
    return fill(f"{first['metric']}: {first['dataset']} AUC", width=56)


def comparison_metadata(auc: pd.DataFrame) -> dict[str, dict[str, str]]:
    data = auc[["metric", "dataset", "context"]].drop_duplicates().copy()
    data["comparison"] = (
        data["metric"].astype(str)
        + ": "
        + data["dataset"].astype(str)
        + " - "
        + data["context"].astype(str)
    )
    return data.set_index("comparison").to_dict("index")


def comparison_datasets(auc: pd.DataFrame) -> dict[str, str]:
    return {
        comparison: str(meta["dataset"])
        for comparison, meta in comparison_metadata(auc).items()
    }


def select_comparisons(
    comparisons: list[str],
    metadata: Mapping[str, Mapping[str, str]],
    dataset: str,
) -> list[str]:
    return sorted(
        [name for name in comparisons if metadata.get(name, {}).get("dataset") == dataset],
        key=lambda name: str(metadata[name]["context"]),
    )


def select_comparison_pairs(
    comparisons: list[str],
    metadata: Mapping[str, Mapping[str, str]],
    requested_pairs: Iterable[Mapping[str, str]] | None,
) -> list[tuple[str, str]]:
    if not requested_pairs:
        return list(combinations(comparisons, 2))

    pairs: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for requested in requested_pairs:
        x_dataset = str(requested["x_dataset"])
        y_dataset = str(requested["y_dataset"])
        x_names = select_comparisons(comparisons, metadata, x_dataset)
        y_names = select_comparisons(comparisons, metadata, y_dataset)

        for y_name, x_name in product(y_names, x_names):
            pair = (x_name, y_name)
            if pair[0] == pair[1] or pair in seen:
                continue
            pairs.append(pair)
            seen.add(pair)

    return pairs


def plot_auc_pairs(
    auc: pd.DataFrame,
    output_pdf: Path,
    output_png: Path | None = None,
    color_map: dict[str, str] | None = None,
    comparison_pairs: Iterable[Mapping[str, str]] | None = None,
    n_cols: int = 3,
    subplot_width: float = 5.6,
    subplot_height: float = 4.9,
    font_scale: float = 1.55,
    point_size: float = 96,
    log_scale: bool = True,
) -> None:
    matrix = make_comparison_matrix(auc, value_col="value" if "value" in auc else "auc")
    comparisons = list(matrix.columns)
    metadata = comparison_metadata(auc)
    requested_pairs = list(comparison_pairs or [])
    pairs = select_comparison_pairs(comparisons, metadata, requested_pairs)
    if not pairs:
        raise ValueError("Need at least two contexts to make AUC scatter plots")

    colors = {**DEFAULT_COLORS, **(color_map or {})}
    use_context_grid = len(requested_pairs) == 1
    if use_context_grid:
        x_names = select_comparisons(comparisons, metadata, str(requested_pairs[0]["x_dataset"]))
        y_names = select_comparisons(comparisons, metadata, str(requested_pairs[0]["y_dataset"]))
        n_cols = len(x_names)
        n_rows = len(y_names)
    else:
        n_cols = min(n_cols, len(pairs))
        n_rows = ceil(len(pairs) / n_cols)

    sns.set_theme(context="notebook", style="whitegrid", font_scale=font_scale)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(subplot_width * n_cols, subplot_height * n_rows),
        squeeze=False,
        constrained_layout=True,
        sharex=use_context_grid,
        sharey=use_context_grid,
    )

    for idx, (ax, (x_name, y_name)) in enumerate(zip(axes.flat, pairs)):
        row_idx = idx // n_cols
        col_idx = idx % n_cols
        plot_frame = matrix[[x_name, y_name]].dropna().reset_index()
        if log_scale:
            plot_frame = plot_frame[
                plot_frame[x_name].gt(0) & plot_frame[y_name].gt(0)
            ].copy()
        point_colors = [colors.get(lnp, "#333333") for lnp in plot_frame["lnp"]]

        if len(plot_frame) >= 2:
            if log_scale:
                slope, intercept = np.polyfit(
                    np.log10(plot_frame[x_name]),
                    np.log10(plot_frame[y_name]),
                    deg=1,
                )
                x_fit = np.geomspace(plot_frame[x_name].min(), plot_frame[x_name].max(), 200)
                y_fit = 10 ** (intercept + slope * np.log10(x_fit))
                ax.plot(x_fit, y_fit, color="#111111", linewidth=1.9, zorder=2)
            else:
                sns.regplot(
                    data=plot_frame,
                    x=x_name,
                    y=y_name,
                    scatter=False,
                    color="#111111",
                    line_kws={"linewidth": 1.9},
                    ci=95,
                    ax=ax,
                )

        ax.scatter(
            plot_frame[x_name],
            plot_frame[y_name],
            c=point_colors,
            s=point_size,
            linewidth=0.85,
            edgecolor="white",
            zorder=3,
        )
        if log_scale:
            ax.set_xscale("log")
            ax.set_yscale("log")
            for axis, subs in ((ax.xaxis, (1.0, 5.0)), (ax.yaxis, (1.0, 2.0, 5.0))):
                axis.set_major_locator(LogLocator(base=10, subs=subs))
                axis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
                axis.set_minor_locator(NullLocator())
            if metadata[y_name]["dataset"] == "in vitro protein expression":
                ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))

        rho, p_value, n = spearman_pair(plot_frame[x_name], plot_frame[y_name])
        annotation = f"Spearman rho = {rho:.2f}\np = {p_value:.2g}; n = {n}"
        if np.isnan(rho):
            annotation = f"Spearman rho = NA\nn = {n}"
        ax.text(
            0.04,
            0.96,
            annotation,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=12,
            bbox={"boxstyle": "round,pad=0.32", "facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
        )
        if use_context_grid:
            ax.set_xlabel("")
            ax.set_ylabel("")
            if row_idx == 0:
                ax.set_title(str(metadata[x_name]["context"]), fontsize=16, pad=12)
            if col_idx == n_cols - 1:
                ax.annotate(
                    str(metadata[y_name]["context"]),
                    xy=(1.04, 0.5),
                    xycoords="axes fraction",
                    rotation=270,
                    ha="left",
                    va="center",
                    fontsize=16,
                )
        else:
            ax.set_xlabel(comparison_label(x_name))
            ax.set_ylabel(comparison_label(y_name))
        ax.tick_params(axis="both", labelsize=12)
        ax.margins(x=0.08, y=0.08)
        ax.grid(True, color="#E5E5E5", linewidth=0.8, which="major")
        ax.grid(True, color="#F1F1F1", linewidth=0.5, which="minor")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    if use_context_grid:
        fig.supxlabel(axis_label(x_names, metadata), fontsize=17)
        fig.supylabel(axis_label(y_names, metadata), fontsize=17)

    for ax in axes.flat[len(pairs):]:
        ax.axis("off")

    handles = [
        Line2D(
            [0],
            [0],
            label=str(lnp),
            marker="o",
            markersize=7,
            linestyle="",
            markerfacecolor=colors.get(lnp, "#333333"),
            markeredgecolor="white",
            markeredgewidth=0.8,
        )
        for lnp in matrix.index
    ]
    if handles:
        fig.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0),
            ncols=min(4, len(handles)),
            fontsize=10,
            frameon=False,
        )

    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, bbox_inches="tight")
    if output_png is not None:
        fig.savefig(output_png, dpi=400, bbox_inches="tight")
    plt.close(fig)
