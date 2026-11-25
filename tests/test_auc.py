import math

import pandas as pd

from lucy_bl251.auc import compute_auc, compute_cyto_abundance, make_comparison_matrix


def test_compute_cyto_abundance_allocates_qpcr_by_cytosolic_read_share():
    data = pd.DataFrame(
        {
            "dataset": ["invivo", "invivo", "invivo"],
            "context": ["Spleen", "Spleen", "Spleen"],
            "lnp": ["A", "A", "B"],
            "time_h": [3, 3, 3],
            "rep": [1, 1, 1],
            "embedded_barcode": ["A1", "A2", "B1"],
            "cytosol_reads": [2, 3, 5],
            "total_reads": [10, 10, 10],
            "qpcr_norm_gl": [100, 100, 100],
        }
    )

    result = compute_cyto_abundance(data).sort_values("lnp")

    assert result["lnp"].tolist() == ["A", "B"]
    assert result["cyto_abund"].tolist() == [50, 50]


def test_compute_auc_uses_time_means_before_trapezoid_auc():
    data = pd.DataFrame(
        {
            "dataset": ["d", "d", "d", "d"],
            "context": ["c", "c", "c", "c"],
            "lnp": ["A", "A", "A", "A"],
            "time_h": [0, 0, 2, 2],
            "cyto_abund": [0, 2, 2, 4],
        }
    )

    result = compute_auc(data, "cyto_abund", "cyto")

    assert result.loc[0, "auc"] == 4


def test_compute_auc_supports_a_non_time_x_column():
    data = pd.DataFrame(
        {
            "dataset": ["d", "d"],
            "context": ["c", "c"],
            "lnp": ["A", "A"],
            "dose_ng": [1, 3],
            "protein_expression": [10, 20],
        }
    )

    result = compute_auc(data, "protein_expression", "protein", x_col="dose_ng")

    assert result.loc[0, "auc"] == 30


def test_make_comparison_matrix_pivots_auc_by_lnp():
    auc = pd.DataFrame(
        {
            "metric": ["cyto", "cyto"],
            "dataset": ["in vivo", "in vitro"],
            "context": ["Spleen", "Hep3B"],
            "lnp": ["A", "A"],
            "auc": [1.5, 2.5],
        }
    )

    matrix = make_comparison_matrix(auc)

    assert matrix.loc["A", "cyto: in vivo - Spleen"] == 1.5
    assert matrix.loc["A", "cyto: in vitro - Hep3B"] == 2.5
    assert not math.isnan(matrix.loc["A", "cyto: in vivo - Spleen"])
