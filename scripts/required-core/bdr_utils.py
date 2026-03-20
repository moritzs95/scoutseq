#!/usr/bin/env python3

"""
Shared helpers for BD Rhapsody barcode and sample-tag normalization.
"""

from __future__ import annotations

import csv
import os
import re
from functools import lru_cache
from pathlib import Path

import pandas as pd


BDR_BARCODE_PATTERN = re.compile(
    r"^(?:A|GT|TCA)?(?P<cls1>[ACGT]{9})GTGA(?P<cls2>[ACGT]{9})GACA(?P<cls3>[ACGT]{9})$"
)


def resolve_bdr_resource_dir(resource_dir: str | Path | None = None) -> Path:
    """Resolve the folder that contains the BDR integration assets."""
    if resource_dir:
        return Path(resource_dir).expanduser().resolve()

    env_dir = os.environ.get("SCOUT_BDR_DIR", "").strip()
    if env_dir:
        return Path(env_dir).expanduser().resolve()

    return Path(__file__).resolve().parents[2] / "BDR_integration"


def hamming_distance(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Barcode segments must have equal length.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(left, right))


def label_to_index96(label: str) -> str:
    cls1, cls2, cls3 = (int(value) for value in label.split("-"))
    return str((cls1 - 1) * 96 * 96 + (cls2 - 1) * 96 + cls3)


@lru_cache(maxsize=None)
def _load_cls_lookup(resource_dir: str) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    base_dir = Path(resource_dir)

    def read_lookup(filename: str) -> dict[str, int]:
        with (base_dir / filename).open() as handle:
            return {
                line.strip(): index
                for index, line in enumerate(handle, start=1)
                if line.strip()
            }

    return (
        read_lookup("BD_CLS1.txt"),
        read_lookup("BD_CLS2.txt"),
        read_lookup("BD_CLS3.txt"),
    )


def _resolve_segment_index(segment: str, lookup: dict[str, int], max_distance: int = 1) -> int | None:
    if segment in lookup:
        return lookup[segment]

    for known_barcode, index in lookup.items():
        if hamming_distance(segment, known_barcode) <= max_distance:
            return index

    return None


def normalize_bdr_barcode(raw_barcode: object, resource_dir: str | Path | None = None) -> str:
    """Convert a raw BD Rhapsody CLS barcode into its numeric Cell_Index."""
    if pd.isna(raw_barcode):
        return ""

    barcode = str(raw_barcode).strip()
    if not barcode:
        return ""
    if barcode.isdigit():
        return barcode

    match = BDR_BARCODE_PATTERN.match(barcode)
    if not match:
        return barcode

    cls1_lookup, cls2_lookup, cls3_lookup = _load_cls_lookup(str(resolve_bdr_resource_dir(resource_dir)))
    cls1_index = _resolve_segment_index(match.group("cls1"), cls1_lookup)
    cls2_index = _resolve_segment_index(match.group("cls2"), cls2_lookup)
    cls3_index = _resolve_segment_index(match.group("cls3"), cls3_lookup)

    if None in (cls1_index, cls2_index, cls3_index):
        return barcode

    return label_to_index96(f"{cls1_index}-{cls2_index}-{cls3_index}")


def normalize_bdr_barcode_series(series: pd.Series, resource_dir: str | Path | None = None) -> pd.Series:
    return series.apply(lambda value: normalize_bdr_barcode(value, resource_dir=resource_dir))


def load_bdr_sample_tags(resource_dir: str | Path | None = None) -> pd.DataFrame:
    """Read BD sample-tag calls while skipping the comment preamble lines."""
    sample_tag_path = resolve_bdr_resource_dir(resource_dir) / "BDR_Sample_Tag_Calls.csv"
    if not sample_tag_path.exists():
        return pd.DataFrame(columns=["bc", "Sample_Tag", "Sample_Name"])

    rows: list[dict[str, str]] = []
    with sample_tag_path.open(newline="") as handle:
        reader = csv.DictReader(line for line in handle if not line.startswith("#"))
        for row in reader:
            rows.append(
                {
                    "bc": str(row["Cell_Index"]).strip(),
                    "Sample_Tag": str(row["Sample_Tag"]).strip(),
                    "Sample_Name": str(row["Sample_Name"]).strip(),
                }
            )

    return pd.DataFrame(rows)
