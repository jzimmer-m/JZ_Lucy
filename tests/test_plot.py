from lucy_bl251.plot import select_comparison_pairs


def test_select_comparison_pairs_uses_requested_dataset_pairs():
    comparisons = ["cyto: in vivo - Liver", "cyto: in vitro - Hep3B", "protein: in vitro protein - Jurkat"]
    metadata = {
        "cyto: in vivo - Liver": {"dataset": "in vivo", "context": "Liver"},
        "cyto: in vitro - Hep3B": {"dataset": "in vitro", "context": "Hep3B"},
        "protein: in vitro protein - Jurkat": {"dataset": "protein", "context": "Jurkat"},
    }

    result = select_comparison_pairs(
        comparisons,
        metadata,
        [{"x_dataset": "in vitro", "y_dataset": "protein"}],
    )

    assert result == [("cyto: in vitro - Hep3B", "protein: in vitro protein - Jurkat")]


def test_select_comparison_pairs_arranges_y_contexts_as_rows_and_x_contexts_as_columns():
    comparisons = [
        "cyto: in vitro - B",
        "cyto: in vitro - A",
        "cyto: in vivo - Y",
        "cyto: in vivo - X",
    ]
    metadata = {
        "cyto: in vitro - A": {"dataset": "in vitro", "context": "A"},
        "cyto: in vitro - B": {"dataset": "in vitro", "context": "B"},
        "cyto: in vivo - X": {"dataset": "in vivo", "context": "X"},
        "cyto: in vivo - Y": {"dataset": "in vivo", "context": "Y"},
    }

    result = select_comparison_pairs(
        comparisons,
        metadata,
        [{"x_dataset": "in vitro", "y_dataset": "in vivo"}],
    )

    assert result == [
        ("cyto: in vitro - A", "cyto: in vivo - X"),
        ("cyto: in vitro - B", "cyto: in vivo - X"),
        ("cyto: in vitro - A", "cyto: in vivo - Y"),
        ("cyto: in vitro - B", "cyto: in vivo - Y"),
    ]
