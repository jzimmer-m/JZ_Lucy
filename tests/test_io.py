import pandas as pd

from lucy_bl251.io import apply_lnp_aliases, build_cyto_input


def test_build_cyto_input_marks_repeated_lnp_plexes_for_downstream_filtering(tmp_path):
    scar_df = tmp_path / "scar_df.csv"
    rna_metadata = tmp_path / "rna_metadata.csv"
    sample_metadata = tmp_path / "sample_metadata.csv"
    qpcr = tmp_path / "qpcr.csv"

    pd.DataFrame(
        {
            "flowcell": ["fc", "fc", "fc"],
            "barcode": ["barcode01", "barcode01", "barcode01"],
            "embedded_barcode": ["BC1", "BC2", "BC3"],
            "Cytosol": [10, 20, 30],
            "Endosome": [1, 2, 3],
            "total": [11, 22, 33],
        }
    ).to_csv(scar_df, index=False)
    pd.DataFrame(
        {
            "internal_barcode": ["BC1", "BC2", "BC3"],
            "LNP": ["86/56", "86/56", "716/56"],
            "plex": [20, 125, 20],
        }
    ).to_csv(rna_metadata, index=False)
    pd.DataFrame(
        {
            "sample": ["sample-a"],
            "barcode": ["barcode01"],
            "flowcell": ["fc"],
            "tissue": ["Spleen"],
            "time": [3],
            "rep": [1],
        }
    ).to_csv(sample_metadata, index=False)
    pd.DataFrame({"sample": ["sample-a"], "norm_GL": [100]}).to_csv(qpcr, index=False)

    result = build_cyto_input(
        scar_df_path=scar_df,
        rna_metadata_path=rna_metadata,
        sample_metadata_path=sample_metadata,
        qpcr_path=qpcr,
        dataset="In vivo",
        context_col="tissue",
        sample_join_cols=["sample"],
        keep_lnp_plexes={"86/56": [20]},
    )

    assert result["lnp"].tolist() == ["86/56 20", "86/56 125", "716/56"]
    assert result["include_lnp"].tolist() == [True, False, True]


def test_apply_lnp_aliases_uses_canonical_names():
    frame = pd.DataFrame({"lnp": [" 86/56 11-plex ", "716/56"], "value": [1, 2]})

    result = apply_lnp_aliases(frame, {"86/56 11-plex": "86/56"})

    assert result["lnp"].tolist() == ["86/56", "716/56"]
