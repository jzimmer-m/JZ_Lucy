import json

from lucy_bl251.lucy_studio import normalize_lucy_studio_imaging_export


def test_normalize_lucy_studio_imaging_export_joins_layout_to_raw_text(tmp_path):
    raw_text = "\r\n".join(
        [
            "Plate Name\tMHPG2G1",
            "",
            "[Data]",
            "\t".join(
                [
                    "Row",
                    "Column",
                    "Cells Selected - Intensity Cell Alexa 488 Mean - Mean per Well",
                ]
            ),
            "1\t1\t100",
            "1\t2\t999",
            "2\t1\t888",
        ]
    )
    fixture = {
        "plateLayouts": [
            {
                "signature": json.dumps(
                    [
                        "PLAN-IMG-CL1-B1-G1",
                        "compact",
                        "HepG2",
                        0,
                        [
                            ["A1", "FORMULATION", 5, "86/56 11 member pool", "LNP-F13440", 200, 1],
                            ["A2", "CELLS_ONLY", "", ""],
                            [
                                "B1",
                                "FORMULATION",
                                1,
                                "86/56 1000 member pool",
                                "LNP-F13436",
                                200,
                                1,
                            ],
                        ],
                    ]
                )
            }
        ],
        "results": [{"rawText": raw_text}],
    }
    path = tmp_path / "lucy_studio.json"
    path.write_text(json.dumps(fixture), encoding="utf-8")

    result = normalize_lucy_studio_imaging_export(path)

    assert result["context"].tolist() == ["HepG2"]
    assert result["lnp"].tolist() == ["86/56 11 member pool"]
    assert result["dose_ng"].tolist() == [200.0]
    assert result["protein_expression"].tolist() == [100]
