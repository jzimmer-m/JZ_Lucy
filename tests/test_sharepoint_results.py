import pandas as pd

from lucy_bl251.sharepoint_results import normalize_scarseq_result_export


def test_normalize_scarseq_result_export_marks_repeated_plexes_and_excludes_contexts():
    frame = pd.DataFrame(
        {
            "internal_barcode": ["BC1", "BC2", "BC3", "BC4"],
            "Cytosol": [10, 20, 30, 40],
            "total": [11, 22, 33, 44],
            "LNP": ["86/56", "86/56", "716/56", "86/56"],
            "time": [3, 3, 3, 0],
            "rep": [1, 1, 1, 1],
            "norm_GL": [100, 100, 100, None],
            "plex": [20, 125, 20, 20],
            "cell": ["Hep3B", "Hep3B", "Hep3B", "DS"],
        }
    )

    result = normalize_scarseq_result_export(
        frame,
        dataset="in vitro Scar-seq",
        context_col="cell",
        exclude_contexts={"DS"},
        keep_lnp_plexes={"86/56": [20]},
    )

    assert result["lnp"].tolist() == ["86/56 20", "86/56 125", "716/56"]
    assert result["include_lnp"].tolist() == [True, False, True]
    assert result["context"].tolist() == ["Hep3B", "Hep3B", "Hep3B"]
    assert result["embedded_barcode"].tolist() == ["BC1", "BC2", "BC3"]
    assert result["cytosol_reads"].tolist() == [10, 20, 30]
