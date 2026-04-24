#!/usr/bin/env python3

"""
Remove IncorrectHDR alleles from EditingOutcomesByCellBarcode.xlsx, shift the
remaining alleles left, and recompute per-row totals.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def max_allele_index(columns: list[str]) -> int:
    indices = []
    for column in columns:
        match = re.fullmatch(r"Allele(\d+)", column)
        if match:
            indices.append(int(match.group(1)))
    return max(indices, default=0)


def recalc_ploidy(allele_count: int, total_umicount: int) -> str:
    if allele_count <= 0:
        return ""
    if allele_count == 1:
        return "haploid" if total_umicount == 1 else "diploid"
    if allele_count == 2:
        return "diploid"
    if allele_count == 3:
        return "triploid"
    return "polyploid"


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()

    df = pd.read_excel(input_path)
    max_index = max_allele_index(list(df.columns))
    if max_index == 0:
      output_path.parent.mkdir(parents=True, exist_ok=True)
      df.to_excel(output_path, index=False, engine="openpyxl")
      return

    rows: list[dict[str, object]] = []
    allele_fields = [
        ("Allele", "Allele{idx}"),
        ("rc", "rc_{idx}"),
        ("rc_perc", "rc_perc_{idx}"),
        ("UMIcount", "UMIcount_{idx}"),
        ("Sequence", "Sequence_{idx}"),
        ("HDist_HDRBC", "HDist_HDRBC_{idx}"),
        ("AlnRefSeq", "AlnRefSeq_{idx}"),
        ("RepairOutcome", "RepairOutcome_{idx}"),
    ]

    for _, row in df.iterrows():
        kept = []
        for idx in range(1, max_index + 1):
            repair_outcome = row.get(f"RepairOutcome_{idx}")
            allele = row.get(f"Allele{idx}")
            if pd.isna(allele):
                continue
            if str(repair_outcome) == "IncorrectHDR":
                continue
            entry = {name: row.get(template.format(idx=idx)) for name, template in allele_fields}
            kept.append(entry)

        if not kept:
            continue

        updated = row.to_dict()
        total_rc = sum(float(entry["rc"]) for entry in kept if pd.notna(entry["rc"]))
        total_umicount = sum(float(entry["UMIcount"]) for entry in kept if pd.notna(entry["UMIcount"]))
        updated["allele_count"] = len(kept)
        updated["total_rc"] = total_rc
        updated["total_UMIcount"] = total_umicount
        updated["Ploidy"] = recalc_ploidy(len(kept), int(total_umicount))
        updated["DetailedRepairOutcome"] = "_".join(str(entry["Allele"]) for entry in kept)

        for idx in range(1, max_index + 1):
            if idx <= len(kept):
                entry = kept[idx - 1]
                updated[f"Allele{idx}"] = entry["Allele"]
                updated[f"rc_{idx}"] = entry["rc"]
                updated[f"rc_perc_{idx}"] = (float(entry["rc"]) / total_rc * 100) if total_rc else 0
                updated[f"UMIcount_{idx}"] = entry["UMIcount"]
                updated[f"Sequence_{idx}"] = entry["Sequence"]
                updated[f"HDist_HDRBC_{idx}"] = entry["HDist_HDRBC"]
                updated[f"AlnRefSeq_{idx}"] = entry["AlnRefSeq"]
                updated[f"RepairOutcome_{idx}"] = entry["RepairOutcome"]
            else:
                updated[f"Allele{idx}"] = None
                updated[f"rc_{idx}"] = None
                updated[f"rc_perc_{idx}"] = None
                updated[f"UMIcount_{idx}"] = None
                updated[f"Sequence_{idx}"] = None
                updated[f"HDist_HDRBC_{idx}"] = None
                updated[f"AlnRefSeq_{idx}"] = None
                updated[f"RepairOutcome_{idx}"] = None

        rows.append(updated)

    output_df = pd.DataFrame(rows, columns=df.columns)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_excel(output_path, index=False, engine="openpyxl")


if __name__ == "__main__":
    main()
