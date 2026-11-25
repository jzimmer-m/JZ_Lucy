from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_EXCLUDED_LNPS = frozenset(
    {
        "86/56 1000 member pool",
        "86/56 500 member pool",
        "86/56 250 member pool",
        "86/56 125 member pool",
    }
)
DEFAULT_CELL_ALIASES = {"Jurkats": "Jurkat"}
DEFAULT_INTENSITY_COL = "Cells Selected - Intensity Cell Alexa 488 Mean - Mean per Well"

PLATE_PATTERNS = (
    (re.compile(r"^MHPG2G(?P<group>\d+)$"), "HepG2"),
    (re.compile(r"^MSKH1G(?P<group>\d+)$"), "SK-Hep-1"),
    (re.compile(r"^MJRKTG(?P<group>\d+)$"), "Jurkats"),
    (re.compile(r"^MTHPMG(?P<group>\d+)$"), "THP-1"),
)

PROTEIN_COLUMNS = [
    "dataset",
    "context",
    "lnp",
    "dose_ng",
    "rep",
    "protein_expression",
    "plate",
    "well",
    "fm_id",
]


@dataclass(frozen=True)
class WellLayout:
    lnp: str
    fm_id: str
    dose_ng: float
    rep: int


def iter_property_values(node: Any, name: str) -> Iterable[Any]:
    if isinstance(node, dict):
        for key, value in node.items():
            if key == name:
                yield value
            yield from iter_property_values(value, name)
    elif isinstance(node, list):
        for value in node:
            yield from iter_property_values(value, name)


def parse_plate_name(plate_name: str) -> tuple[str, str]:
    for pattern, cell in PLATE_PATTERNS:
        match = pattern.match(plate_name)
        if match:
            return cell, f"G{match.group('group')}"
    raise ValueError(f"Unrecognized LucyStudio plate name: {plate_name}")


def extract_plate_name(raw_text: str) -> str:
    match = re.search(r"^Plate Name\t(?P<plate>[^\r\n]+)", raw_text, flags=re.MULTILINE)
    if not match:
        raise ValueError("LucyStudio raw text is missing a Plate Name line")
    return match.group("plate")


def parse_raw_text_table(raw_text: str) -> pd.DataFrame:
    sections = re.split(r"\[Data\]\r?\n", raw_text, maxsplit=1)
    if len(sections) != 2:
        raise ValueError("LucyStudio raw text is missing a [Data] section")
    return pd.read_csv(StringIO(sections[1]), sep="\t")


def well_label(row: pd.Series) -> str:
    return f"{chr(64 + int(row['Row']))}{int(row['Column'])}"


def extract_layouts(document: Mapping[str, Any]) -> dict[tuple[str, str], dict[str, WellLayout]]:
    layouts = {}
    for signature_text in iter_property_values(document, "signature"):
        signature = json.loads(signature_text)
        if len(signature) < 5:
            continue

        plan_name = str(signature[0])
        group_match = re.search(r"-(G\d+)$", plan_name)
        if not group_match:
            continue

        cell = str(signature[2])
        wells = {}
        for entry in signature[4]:
            if len(entry) < 7 or entry[1] != "FORMULATION":
                continue
            wells[str(entry[0])] = WellLayout(
                lnp=str(entry[3]).strip(),
                fm_id=str(entry[4]).strip(),
                dose_ng=float(entry[5]),
                rep=int(entry[6]),
            )
        layouts[(cell, group_match.group(1))] = wells
    return layouts


def normalize_lucy_studio_imaging_export(
    path: Path,
    *,
    dataset: str = "in vitro protein expression",
    value_col: str = DEFAULT_INTENSITY_COL,
    exclude_lnps: Iterable[str] = DEFAULT_EXCLUDED_LNPS,
    cell_aliases: Mapping[str, str] = DEFAULT_CELL_ALIASES,
) -> pd.DataFrame:
    with path.open("r", encoding="utf-8") as handle:
        document = json.load(handle)

    layouts = extract_layouts(document)
    excluded = {lnp.strip() for lnp in exclude_lnps}
    rows = []
    for raw_text in iter_property_values(document, "rawText"):
        plate = extract_plate_name(raw_text)
        cell, group = parse_plate_name(plate)
        layout = layouts.get((cell, group))
        if layout is None:
            raise ValueError(f"No LucyStudio layout found for {cell} {group}")

        data = parse_raw_text_table(raw_text)
        if value_col not in data:
            raise ValueError(f"{plate} is missing the {value_col!r} measurement column")

        for _, row in data.iterrows():
            if pd.isna(row["Row"]) or pd.isna(row["Column"]):
                continue

            well = well_label(row)
            well_layout = layout.get(well)
            if well_layout is None or not well_layout.lnp or well_layout.lnp in excluded:
                continue

            rows.append(
                {
                    "dataset": dataset,
                    "context": cell_aliases.get(cell, cell),
                    "lnp": well_layout.lnp,
                    "dose_ng": well_layout.dose_ng,
                    "rep": well_layout.rep,
                    "protein_expression": row[value_col],
                    "plate": plate,
                    "well": well,
                    "fm_id": well_layout.fm_id,
                }
            )

    if not rows:
        raise ValueError(f"No protein-expression rows were parsed from {path}")

    result = pd.DataFrame(rows, columns=PROTEIN_COLUMNS)
    result["dose_ng"] = pd.to_numeric(result["dose_ng"], errors="coerce")
    result["rep"] = pd.to_numeric(result["rep"], errors="coerce")
    result["protein_expression"] = pd.to_numeric(result["protein_expression"], errors="coerce")
    return result.sort_values(["context", "lnp", "dose_ng", "rep", "plate", "well"]).reset_index(
        drop=True
    )
